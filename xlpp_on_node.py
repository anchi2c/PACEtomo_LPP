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

def _box_std(image, box_size=None, box_center=None):
    """Standard deviation of a square box of the image. Default box is centered half of the image."""
    image = np.asarray(image, dtype=float)
    if box_size is None:
        box_size = min(image.shape) // 2
    if box_center is None:
        box_center = (image.shape[0] // 2, image.shape[1] // 2)
    half = box_size // 2
    r0 = max(0, int(box_center[0]) - half)
    c0 = max(0, int(box_center[1]) - half)
    return float(image[r0:r0 + box_size, c0:c0 + box_size].std())

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

def find_maximum_contrast_xts(number_of_changes, box_size=None, box_center=None, settle_delay_s=0.0,
                              frequency=1.0, show_plot=True):
    """
    Scan XLensDeflector(2) on a number_of_changes x number_of_changes grid
    oriented along the two row vectors of ronchiCorrMatrix (xt change per rad
    of phase of each k vector, as solved by ronchi_cal), covering one
    phase period (2 pi) along each. Acquire a Trial image at each point, take
    the standard deviation of a boxed area as the contrast, and fit a 2D
    cosine to find the x-tilt of maximum contrast.
    frequency: expected contrast periods per 2 pi of phase, used to start the fit.
    show_plot: show std_grid heatmaps with the fit (display_util.showContrastGrid).
    XLensDeflector(2) is restored to its starting value at the end.
    Returns (best_xt, fit_params, std_grid, xt_grid).
    """
    corr = np.asarray(ronchiCorrMatrix, dtype=float) * corr_scale
    xt0 = np.array(ronchi_cal.sem.ReportXLensDeflector(2), dtype=float)
    # phase coordinates centered on the start, one period per axis
    phases = 2 * np.pi * (np.arange(number_of_changes) - (number_of_changes - 1) / 2.0) / number_of_changes
    std_grid = np.zeros((number_of_changes, number_of_changes))
    xt_grid = np.zeros((number_of_changes, number_of_changes, 2))
    try:
        for i, u in enumerate(phases):
            for j, v in enumerate(phases):
                # row vector of ronchiCorrMatrix corresponds to one k vector
                xt = xt0 + u * corr[0] + v * corr[1]
                xt_grid[i, j] = xt
                ronchi_cal.sem.SetXLensDeflector(2, xt[0], xt[1])
                if settle_delay_s > 0:
                    ronchi_cal.sem.Delay(settle_delay_s, "s")
                ronchi_cal.sem.T()
                image = np.asarray(ronchi_cal.sem.bufferImage("A"))
                std_grid[i, j] = _box_std(image, box_size, box_center)
                print(f'phase ({u:.2f}, {v:.2f}) rad xt ({xt[0]:.6e}, {xt[1]:.6e}) std {std_grid[i, j]:.3f}')
    finally:
        ronchi_cal.sem.SetXLensDeflector(2, xt0[0], xt0[1])
    u_grid, v_grid = np.meshgrid(phases, phases, indexing='ij')
    params, u_max, v_max = fit_cosine_2d(u_grid, v_grid, std_grid, frequency)
    best_xt = xt0 + u_max * corr[0] + v_max * corr[1]
    c, a_u, f_u, p_u, a_v, f_v, p_v = params
    print(f'cosine fit: offset {c:.3f}, amplitude u {a_u:.3f} v {a_v:.3f}, frequency u {f_u:.3f} v {f_v:.3f}')
    print(f'maximum contrast at phase ({u_max:.3f}, {v_max:.3f}) rad from start')
    print(f'maximum contrast xt: ({best_xt[0]:.6e}, {best_xt[1]:.6e}); start xt ({xt0[0]:.6e}, {xt0[1]:.6e})')
    display_util.addContrastGrid(phases, std_grid, params, u_max, v_max,
                                 subtitle=f'{number_of_changes}x{number_of_changes}')
    if show_plot:
        display_util.showContrastGrid()
    return best_xt, params, std_grid, xt_grid

if __name__=='__main__':
    # ronchiStartC3 calibration
    if ronchi_cal.ronchi_sem_lib.is_simu == True:
        # offset from on-plane to simulate constrast change
        ronchi_cal.ronchi_sem_lib.ronchiC3Offset = 1
    best_xt, params, std_grid, xt_grid = find_maximum_contrast_xts(5, box_size=None, box_center=None, settle_delay_s=0.0,
                              frequency=1.0, show_plot=True)
    print(f'Set X-tilt to {best_xt[0]:.6f}, {best_xt[1]:.6f}')
    ronchi_cal.sem.SetXLensDeflector(2, best_xt[0], best_xt[1])
    
