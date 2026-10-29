#!Python
import ronchi_cal
import display_util
import numpy as np

ronchi_cal.ronchi_sem_lib.checkRonchigramSetup()
pixel_size_angs = ronchi_cal.sem.ReportCurrentPixelSize('T')
xt_tilt = 1 #(scaled at 1e-5 rad)
corr_scale = 1e-5
image_binning = 2
ronchi_binning = 32
# The same value has always been used.
ronchi_fixed_pixel_size_um = 0.98e-4 * 2

ronchiCorrMatrix = ronchi_cal.ronchi_sem_lib.ronchiCorrMatrix
# quadrants selected by the last find_maximum_contrast_xts, reused in refinement
last_selected_quadrants = None

def _box_normalized_std(image, box_size=None, box_center=None):
    """
    Normalized standard deviation (std / mean) of a square box of the image.
    Default box is centered half of the image.
    """
    image = np.asarray(image, dtype=float)
    if box_size is None:
        box_size = min(image.shape) // 2
    if box_center is None:
        box_center = (image.shape[0] // 2, image.shape[1] // 2)
    half = box_size // 2
    r0 = max(0, int(box_center[0]) - half)
    c0 = max(0, int(box_center[1]) - half)
    box = image[r0:r0 + box_size, c0:c0 + box_size]
    return float(box.std() / box.mean())

QUADRANT_NAMES = ('top-left', 'top-right', 'bottom-left', 'bottom-right')
QUADRANT_STATS = ('mean', 'min', 'max', 'std')

def _quadrant_stats(image, box_size=None, box_center=None):
    """
    Pixel statistics of the four quadrants of the box used by _box_normalized_std.
    Returns {stat: array of 4} for stat in QUADRANT_STATS plus 'std/mean',
    quadrants in the order of QUADRANT_NAMES (numpy row 0 is top).
    """
    image = np.asarray(image, dtype=float)
    if box_size is None:
        box_size = min(image.shape) // 2
    if box_center is None:
        box_center = (image.shape[0] // 2, image.shape[1] // 2)
    half = box_size // 2
    r0 = max(0, int(box_center[0]) - half)
    c0 = max(0, int(box_center[1]) - half)
    box = image[r0:r0 + box_size, c0:c0 + box_size]
    rm, cm = box.shape[0] // 2, box.shape[1] // 2
    quadrants = (box[:rm, :cm], box[:rm, cm:], box[rm:, :cm], box[rm:, cm:])
    stats = {stat: np.array([getattr(q, stat)() for q in quadrants]) for stat in QUADRANT_STATS}
    stats['std/mean'] = stats['std'] / stats['mean']
    return stats

def _cosine_2d(uv, c, a_u, f_u, p_u, a_v, f_v, p_v):
    u, v = uv
    return c + a_u * np.cos(f_u * u - p_u) + a_v * np.cos(f_v * v - p_v)

def fit_cosine_2d(u, v, values, frequency=1.0):
    """
    Fit values = c + a_u*cos(f_u*u - p_u) + a_v*cos(f_v*v - p_v).
    u, v are phase coordinates (rad) along the two ronchiCorrMatrix vectors.
    A linear least squares fit at f = frequency (1 is one contrast period per
    2 pi of phase) gives the starting values; a nonlinear fit then refines
    frequencies and phases.
    Returns (params, u_max, v_max) where (u_max, v_max) is the maximum nearest
    to the scan center.
    """
    u = np.asarray(u, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    values = np.asarray(values, dtype=float).ravel()
    # cos(f*u - p) = cos(p)cos(f*u) + sin(p)sin(f*u) makes the fixed f model linear
    fu, fv = frequency * u, frequency * v
    design = np.column_stack([np.ones_like(u), np.cos(fu), np.sin(fu), np.cos(fv), np.sin(fv)])
    (c, cu, su, cv, sv), *_ = np.linalg.lstsq(design, values, rcond=None)
    params = [c, np.hypot(cu, su), frequency, np.arctan2(su, cu), np.hypot(cv, sv), frequency, np.arctan2(sv, cv)]
    try:
        from scipy.optimize import curve_fit
        refined, _ = curve_fit(_cosine_2d, (u, v), values, p0=params, maxfev=10000)
        # accept only if frequencies stay reasonable
        if 0.25 < abs(refined[2]) < 4 and 0.25 < abs(refined[5]) < 4:
            params = list(refined)
        else:
            print('WARNING: nonlinear cosine fit frequencies unreasonable, using linear fit')
    except Exception as e:
        print(f'WARNING: nonlinear cosine fit failed, using linear fit: {e}')
    c, a_u, f_u, p_u, a_v, f_v, p_v = params
    # make amplitudes and frequencies positive so the maximum is at f*u - p = 0
    if f_u < 0:
        f_u, p_u = -f_u, -p_u
    if a_u < 0:
        a_u, p_u = -a_u, p_u + np.pi
    if f_v < 0:
        f_v, p_v = -f_v, -p_v
    if a_v < 0:
        a_v, p_v = -a_v, p_v + np.pi
    # maximum nearest to the scan center (u = v = 0)
    u_max = (np.mod(p_u + np.pi, 2 * np.pi) - np.pi) / f_u
    v_max = (np.mod(p_v + np.pi, 2 * np.pi) - np.pi) / f_v
    return (c, a_u, f_u, p_u, a_v, f_v, p_v), u_max, v_max

def _scan_contrast_grid(center_xt, phases, corr, box_size=None, box_center=None, settle_delay_s=0.0):
    """
    Acquire a Trial image at each xt = center_xt + u*corr[0] + v*corr[1] for u, v
    in phases and return (contrast_grid, xt_grid, quadrant_stats) indexed [u, v].
    contrast_grid is the normalized std (std / mean) of the whole box;
    quadrant_stats[stat][u, v] has the stat of each box quadrant in the order
    of QUADRANT_NAMES.
    XLensDeflector(2) is restored to its value before the scan.
    """
    xt_before = np.array(ronchi_cal.sem.ReportXLensDeflector(2), dtype=float)
    n = len(phases)
    contrast_grid = np.zeros((n, n))
    quadrant_stats = {stat: np.zeros((n, n, 4)) for stat in QUADRANT_STATS + ('std/mean',)}
    xt_grid = np.zeros((n, n, 2))
    try:
        for i, u in enumerate(phases):
            for j, v in enumerate(phases):
                # row vector of ronchiCorrMatrix corresponds to one k vector
                xt = center_xt + u * corr[0] + v * corr[1]
                xt_grid[i, j] = xt
                ronchi_cal.sem.SetXLensDeflector(2, xt[0], xt[1])
                if settle_delay_s > 0:
                    ronchi_cal.sem.Delay(settle_delay_s, "s")
                ronchi_cal.sem.T()
                image = np.asarray(ronchi_cal.sem.bufferImage("A"))
                contrast_grid[i, j] = _box_normalized_std(image, box_size, box_center)
                for stat, values in _quadrant_stats(image, box_size, box_center).items():
                    quadrant_stats[stat][i, j] = values
                print(f'phase ({u:.2f}, {v:.2f}) rad xt ({xt[0]:.6e}, {xt[1]:.6e}) std/mean {contrast_grid[i, j]:.4f}')
    finally:
        ronchi_cal.sem.SetXLensDeflector(2, xt_before[0], xt_before[1])
    return contrast_grid, xt_grid, quadrant_stats

def _select_low_std_quadrants(quadrant_std_grid, number_of_quadrants=2):
    """
    Pick the quadrants with the lowest median of quadrant_std_grid (std/mean)
    over the whole scan, so quadrants with ice chunks or other strong features
    are excluded, and return (contrast_grid, selected) where contrast_grid is
    the mean of the selected quadrants. The same quadrants are used for every
    grid point.
    """
    medians = np.median(quadrant_std_grid.reshape(-1, 4), axis=0)
    selected = [int(q) for q in np.argsort(medians)[:number_of_quadrants]]
    print('quadrant median std/mean: ' + ', '.join(f'{QUADRANT_NAMES[q]} {medians[q]:.4f}' for q in range(4)))
    print('using quadrants: ' + ', '.join(QUADRANT_NAMES[q] for q in selected))
    return quadrant_std_grid[:, :, selected].mean(axis=2), selected

def fit_quadratic_2d(u, v, values):
    """
    Fit values = a + b*u + c*v + d*u^2 + e*u*v + f*v^2 by least squares.
    Returns (coeffs, u_max, v_max, is_maximum) where (u_max, v_max) is the
    stationary point and is_maximum tells if it is a maximum (negative definite).
    """
    u = np.asarray(u, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    values = np.asarray(values, dtype=float).ravel()
    design = np.column_stack([np.ones_like(u), u, v, u * u, u * v, v * v])
    coeffs, *_ = np.linalg.lstsq(design, values, rcond=None)
    a, b, c, d, e, f = coeffs
    hessian = np.array([[2 * d, e], [e, 2 * f]])
    is_maximum = bool(np.all(np.linalg.eigvalsh(hessian) < 0))
    try:
        u_max, v_max = np.linalg.solve(hessian, [-b, -c])
    except np.linalg.LinAlgError:
        u_max, v_max, is_maximum = 0.0, 0.0, False
    return tuple(coeffs), float(u_max), float(v_max), is_maximum

def find_maximum_contrast_xts(number_of_changes, box_size=None, box_center=None, settle_delay_s=0.0,
                              frequency=1.0, show_plot=True, number_of_quadrants=4):
    """
    Scan x-tilt on a number_of_changes x number_of_changes grid
    oriented along the two row vectors of ronchiCorrMatrix (xt change per rad
    of phase of each k vector), covering one
    phase period (2 pi) along each. Acquire a Trial image at each point, take
    the normalized standard deviation (std / mean) of a boxed area as the
    contrast, and fit a 2D cosine to find the x-tilt of maximum contrast.
    number_of_quadrants: if less than 4, the box is split into four quadrants
    and only the number_of_quadrants with the lowest median std/mean over the
    scan are used as the contrast. Default 4 uses the whole box.
    frequency: expected contrast periods per 2 pi of phase, used to start the fit.
    show_plot: add contrast (std/mean) heatmaps with the fit (display_util.showContrastGrid)
    and quadrant statistics (display_util.showQuadrantStats) for display.
    XLensDeflector(2) is X-tilt and is restored to its starting value at the end.
    Returns (best_xt, fit_params, contrast_grid, xt_grid).
    """
    corr = np.asarray(ronchiCorrMatrix, dtype=float) * corr_scale
    xt0 = np.array(ronchi_cal.sem.ReportXLensDeflector(2), dtype=float)
    # phase coordinates centered on the start, one period per axis
    phases = 2 * np.pi * (np.arange(number_of_changes) - (number_of_changes - 1) / 2.0) / number_of_changes
    box_contrast_grid, xt_grid, quadrant_stats = _scan_contrast_grid(xt0, phases, corr, box_size, box_center, settle_delay_s)
    global last_selected_quadrants
    if number_of_quadrants < 4:
        contrast_grid, last_selected_quadrants = _select_low_std_quadrants(quadrant_stats['std/mean'], number_of_quadrants)
    else:
        contrast_grid = box_contrast_grid
        last_selected_quadrants = None
    u_grid, v_grid = np.meshgrid(phases, phases, indexing='ij')
    params, u_max, v_max = fit_cosine_2d(u_grid, v_grid, contrast_grid, frequency)
    best_xt = xt0 + u_max * corr[0] + v_max * corr[1]
    c, a_u, f_u, p_u, a_v, f_v, p_v = params
    print(f'cosine fit: offset {c:.3f}, amplitude u {a_u:.3f} v {a_v:.3f}, frequency u {f_u:.3f} v {f_v:.3f}')
    print(f'maximum contrast at phase ({u_max:.3f}, {v_max:.3f}) rad from start')
    print(f'maximum contrast xt: ({best_xt[0]:.6e}, {best_xt[1]:.6e}); start xt ({xt0[0]:.6e}, {xt0[1]:.6e})')
    if show_plot:
        display_util.addContrastGrid(phases, contrast_grid, params, u_max, v_max,
                                 subtitle=f'{number_of_changes}x{number_of_changes}')
        display_util.addQuadrantStats(quadrant_stats, last_selected_quadrants, QUADRANT_NAMES,
                                      subtitle=f'{number_of_changes}x{number_of_changes}')
    return best_xt, params, contrast_grid, xt_grid

def refine_maximum_contrast_xts(best_xt, number_of_changes=5, fraction_of_period=1.0 / 3,
                                box_size=None, box_center=None, settle_delay_s=0.0, show_plot=True,
                                quadrants=None, number_of_quadrants=4):
    """
    Refine best_xt from find_maximum_contrast_xts. Scan a
    number_of_changes x number_of_changes grid centered at best_xt along the
    two row vectors of ronchiCorrMatrix, spanning fraction_of_period of one
    phase period (2 pi) of each k vector. Fit the normalized std (std / mean)
    of the boxed image with a 2D second degree polynomial and take its maximum.
    With quadrants, the contrast is the mean std/mean of the box quadrants given by quadrants
    (indices into QUADRANT_NAMES). If quadrants is None, the quadrants selected
    by the last find_maximum_contrast_xts are used; if there are none, the
    number_of_quadrants with the lowest median std/mean in this scan are used
    (default 4 uses the whole box). Use the same box_size and box_center as
    find_maximum_contrast_xts so the quadrants are the same areas.
    If the fit has no maximum, the best measured grid point is used. A maximum
    outside the grid is clipped to the grid edge.
    XLensDeflector(2) is restored to its value before the scan.
    Returns (refined_xt, coeffs, contrast_grid, xt_grid).
    """
    corr = np.asarray(ronchiCorrMatrix, dtype=float) * corr_scale
    best_xt = np.asarray(best_xt, dtype=float)
    half_span = np.pi * fraction_of_period
    phases = np.linspace(-half_span, half_span, number_of_changes)
    box_contrast_grid, xt_grid, quadrant_stats = _scan_contrast_grid(best_xt, phases, corr, box_size, box_center, settle_delay_s)
    if quadrants is None:
        quadrants = last_selected_quadrants
    if quadrants is not None:
        quadrants = [int(q) for q in quadrants]
        print('using quadrants: ' + ', '.join(QUADRANT_NAMES[q] for q in quadrants))
        contrast_grid = quadrant_stats['std/mean'][:, :, quadrants].mean(axis=2)
    elif number_of_quadrants < 4:
        contrast_grid, quadrants = _select_low_std_quadrants(quadrant_stats['std/mean'], number_of_quadrants)
    else:
        contrast_grid = box_contrast_grid
    u_grid, v_grid = np.meshgrid(phases, phases, indexing='ij')
    coeffs, u_max, v_max, is_maximum = fit_quadratic_2d(u_grid, v_grid, contrast_grid)
    if not is_maximum:
        i, j = np.unravel_index(np.argmax(contrast_grid), contrast_grid.shape)
        u_max, v_max = phases[i], phases[j]
        print(f'WARNING: quadratic fit has no maximum, using best measured grid point ({u_max:.3f}, {v_max:.3f}) rad')
    elif abs(u_max) > half_span or abs(v_max) > half_span:
        print(f'WARNING: quadratic fit maximum ({u_max:.3f}, {v_max:.3f}) rad is outside the grid, clipped to the grid edge')
        u_max = float(np.clip(u_max, -half_span, half_span))
        v_max = float(np.clip(v_max, -half_span, half_span))
    refined_xt = best_xt + u_max * corr[0] + v_max * corr[1]
    print(f'quadratic fit coefficients (1, u, v, u^2, uv, v^2): {np.array2string(np.array(coeffs), precision=4)}')
    print(f'refined maximum contrast at phase ({u_max:.3f}, {v_max:.3f}) rad from best_xt')
    print(f'refined maximum contrast xt: ({refined_xt[0]:.6e}, {refined_xt[1]:.6e}); best_xt ({best_xt[0]:.6e}, {best_xt[1]:.6e})')
    if show_plot:
        display_util.addContrastGrid(phases, contrast_grid, coeffs, u_max, v_max,
                                 subtitle=f'refine {number_of_changes}x{number_of_changes}', model='quadratic')
        display_util.addQuadrantStats(quadrant_stats, quadrants, QUADRANT_NAMES,
                                      subtitle=f'refine {number_of_changes}x{number_of_changes}')
    return refined_xt, coeffs, contrast_grid, xt_grid

def phase_distance(xt_a, xt_b):
    """
    Magnitude (rad) of the phase change between two x-tilts, using
    xt = u*corr[0] + v*corr[1] with the row vectors of ronchiCorrMatrix.
    """
    corr = np.asarray(ronchiCorrMatrix, dtype=float) * corr_scale
    delta = np.asarray(xt_b, dtype=float) - np.asarray(xt_a, dtype=float)
    uv = np.linalg.solve(corr.T, delta)
    return float(np.linalg.norm(uv))

if __name__=='__main__':
    # measure xt value at high contrast
    if ronchi_cal.ronchi_sem_lib.is_simu == True:
        # offset from on-plane to simulate constrast change
        ronchi_cal.ronchi_sem_lib.ronchiC3Offset = 1
    best_xt, params, std_grid, xt_grid = find_maximum_contrast_xts(5, box_size=None, box_center=None, settle_delay_s=0.0,
                              frequency=1.0, show_plot=True)
    max_refine_iterations = 5
    # stop refining when best_xt moves less than this phase magnitude (rad)
    refine_phase_threshold = 2 * np.pi / 30
    for i in range(max_refine_iterations):
        fraction = 3+i if i <=1 else 4
        previous_xt = best_xt
        best_xt, coeffs, std_grid, xt_grid = refine_maximum_contrast_xts(best_xt, number_of_changes=5,
                              fraction_of_period=1.0 / (3+i), show_plot=True)
        print(f'Set X-tilt to {best_xt[0]:.6f}, {best_xt[1]:.6f}')
        moved = phase_distance(previous_xt, best_xt)
        print(f'refinement {i + 1}: best_xt moved {moved:.3f} rad (threshold {refine_phase_threshold:.3f} rad)')
        if moved < refine_phase_threshold:
            print(f'refinement converged after {i + 1} iterations')
            break
    else:
        print(f'WARNING: refinement did not converge within {max_refine_iterations} iterations')
    ronchi_cal.sem.SetXLensDeflector(2, best_xt[0], best_xt[1])
    ronchi_cal.sem.T()
    if display_util.Contrast_Phases:
        display_util.showContrastGrid()
    if display_util.Quadrant_Stats:
        display_util.showQuadrantStats()

