#!Python
##############################################################################
# On-cross measurement and correction for xlpp.
#
# The two laser fringe sets cross in the ronchigram. The fringe contrast is
# highest at the center of the cross. Changing the ronchigram C3 offset zooms
# the ronchigram about the beam axis by about 1/|C3 offset|, so a cross that
# is not on the axis moves across a C3 scan. Correcting at each step of a scan
# toward higher zoom brings the cross to the target with increasing precision.
##############################################################################
import numpy as np
import ronchi_cal
import ronchi_lib
import display_util

sem = ronchi_cal.sem
ronchi_sem_lib = ronchi_cal.ronchi_sem_lib

ronchi_sem_lib.checkRonchigramSetup()

def measure_ks_phases(image, pixel_size_um=None, binning=None, peak_radius=None):
    """
    ks (1/um) and phases (rad) of the two fringe sets, the same analysis as
    ronchi_cal.measure_ronchigram_ks_phases on an image already acquired.
    Phases are referenced to the center of the binned image.
    """
    pixel_size_um = ronchi_sem_lib.ronchiPixelSize if pixel_size_um is None else pixel_size_um
    binning = ronchi_sem_lib.ronchiBinning if binning is None else binning
    peak_radius = ronchi_sem_lib.ronchiPeakRadius if peak_radius is None else peak_radius
    binned = ronchi_lib._ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = ronchi_lib._ronchi_find_fourier_centered(binned)
    return ronchi_lib._ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius,
                                             binning=1, fourier_size=image_fft.shape[0])

def _phase_reference(image_shape, binning):
    """Unbinned (row, col) position of the binned image center used as phase reference."""
    binned_shape = np.array(image_shape) // binning
    # binned pixel m covers unbinned m*binning .. m*binning + binning - 1
    return binned_shape / 2.0 * binning + (binning - 1) / 2.0

def fringe_lattice(image_shape, ks, phases, pixel_size_um=None, binning=None):
    """
    Positions (row, col) in unbinned pixels where both fringe phases are a
    multiple of pi: nodes where both are 0 (mod 2 pi) and antinodes where both
    are pi (mod 2 pi), for the fringe model cos(2 pi k.r + phase).
    Returns (positions, lattice_index (a, b), is_node) for points in the image.
    """
    pixel_size_um = ronchi_sem_lib.ronchiPixelSize if pixel_size_um is None else pixel_size_um
    binning = ronchi_sem_lib.ronchiBinning if binning is None else binning
    q = np.asarray(ks, dtype=float) * pixel_size_um   # cycles per unbinned pixel, rows are the two k
    center = _phase_reference(image_shape, binning)
    rows, cols = image_shape
    corners = np.array([[0, 0], [0, cols - 1], [rows - 1, 0], [rows - 1, cols - 1]], dtype=float)
    # phase / pi at the image corners gives the lattice index range
    corner_index = (2 * (corners - center) @ q.T) + np.asarray(phases) / np.pi
    a_range = range(int(np.floor(corner_index[:, 0].min())), int(np.ceil(corner_index[:, 0].max())) + 1)
    b_range = range(int(np.floor(corner_index[:, 1].min())), int(np.ceil(corner_index[:, 1].max())) + 1)
    q_inv = np.linalg.inv(q)
    positions, indices, is_node = [], [], []
    for a in a_range:
        for b in b_range:
            if (a - b) % 2:
                continue    # one phase 0 and the other pi: neither node nor antinode
            r = center + q_inv @ ((np.array([a, b]) - np.asarray(phases) / np.pi) / 2.0)
            if 0 <= r[0] < rows and 0 <= r[1] < cols:
                positions.append(r)
                indices.append((a, b))
                is_node.append(a % 2 == 0)
    return np.array(positions), np.array(indices), np.array(is_node)

def lattice_contrast(image, positions, indices, is_node, box_size, beam_box_size, beam_fraction=0.25):
    """
    Mean of a box_size box centered on each lattice point, and the local fringe
    contrast at each point: node mean - mean of its diagonal antinode
    neighbors, or mean of the neighbor nodes - antinode mean. This cancels slow
    background changes across the ronchigram. Points whose box is outside the
    image, or outside the beam, or without a valid neighbor get nan contrast.
    Outside the beam is where the mean over beam_box_size (about one fringe
    period, which averages out the fringes) is below beam_fraction of the
    bright image level.
    Returns (box_means, contrast).
    """
    image = np.asarray(image, dtype=float)
    rows, cols = image.shape
    half = box_size // 2
    beam_half = beam_box_size // 2
    beam_level = beam_fraction * np.percentile(image, 95)
    box_means = np.full(len(positions), np.nan)
    for n, (r, c) in enumerate(np.round(positions).astype(int)):
        if half <= r < rows - half and half <= c < cols - half:
            beam_mean = image[max(0, r - beam_half):r + beam_half + 1, max(0, c - beam_half):c + beam_half + 1].mean()
            if beam_mean >= beam_level:
                box_means[n] = image[r - half:r + half + 1, c - half:c + half + 1].mean()
    lookup = {tuple(idx): n for n, idx in enumerate(indices)}
    contrast = np.full(len(positions), np.nan)
    for n, (a, b) in enumerate(indices):
        if np.isnan(box_means[n]):
            continue
        neighbors = [lookup.get((a + da, b + db)) for da in (-1, 1) for db in (-1, 1)]
        neighbor_means = [box_means[k] for k in neighbors if k is not None and not np.isnan(box_means[k])]
        if not neighbor_means:
            continue
        sign = 1.0 if is_node[n] else -1.0
        contrast[n] = sign * (box_means[n] - np.mean(neighbor_means))
    return box_means, contrast

def fit_quadratic_2d(x, y, values):
    """
    Fit values = a + b*x + c*y + d*x^2 + e*x*y + f*y^2 by least squares.
    Returns (coeffs, x_max, y_max, is_maximum) where (x_max, y_max) is the
    stationary point and is_maximum tells if it is a maximum.
    """
    x, y, values = (np.asarray(v, dtype=float).ravel() for v in (x, y, values))
    design = np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y])
    coeffs, *_ = np.linalg.lstsq(design, values, rcond=None)
    a, b, c, d, e, f = coeffs
    hessian = np.array([[2 * d, e], [e, 2 * f]])
    is_maximum = bool(np.all(np.linalg.eigvalsh(hessian) < 0))
    try:
        x_max, y_max = np.linalg.solve(hessian, [-b, -c])
    except np.linalg.LinAlgError:
        x_max, y_max, is_maximum = np.nan, np.nan, False
    return tuple(coeffs), float(x_max), float(y_max), is_maximum

def find_cross_center(image, ks, phases, box_fraction=0.25, fit_fraction=0.3, beam_fraction=0.25):
    """
    Center (row, col) of the laser fringe cross in unbinned pixels.
    box_fraction: box size as a fraction of the shorter fringe period.
    fit_fraction: only points with contrast above fit_fraction of the 95th
    percentile contrast are fitted, so the quadratic describes the cross and
    not the flat area outside it.
    Returns a dict with the center, fit and lattice data.
    """
    image = np.asarray(image, dtype=float)
    q = np.asarray(ks, dtype=float) * ronchi_sem_lib.ronchiPixelSize
    periods = 1.0 / np.linalg.norm(q, axis=1)
    box_size = max(3, int(box_fraction * periods.min()))
    positions, indices, is_node = fringe_lattice(image.shape, ks, phases)
    box_means, contrast = lattice_contrast(image, positions, indices, is_node, box_size,
                                           int(periods.max()), beam_fraction)
    valid = ~np.isnan(contrast)
    if valid.sum() < 6:
        raise ValueError(f'only {valid.sum()} valid lattice points, can not fit the cross')
    threshold = fit_fraction * np.percentile(contrast[valid], 95)
    used = valid & (contrast >= threshold)
    if used.sum() < 6:
        used = valid
    coeffs, row_max, col_max, is_maximum = fit_quadratic_2d(positions[used, 0], positions[used, 1], contrast[used])
    if not is_maximum:
        print('WARNING: quadratic fit has no maximum, using the highest contrast point')
        row_max, col_max = positions[valid][np.argmax(contrast[valid])]
    # how much the contrast varies across the image; near 0 when the whole
    # image is inside the cross and its center can not be located
    p5, p95 = np.percentile(contrast[valid], [5, 95])
    relief = float(1 - p5 / p95) if p95 > 0 else 0.0
    inside = bool(0 <= row_max < image.shape[0] and 0 <= col_max < image.shape[1])
    return {
        'center': np.array([row_max, col_max]), 'is_maximum': is_maximum, 'coeffs': coeffs,
        'relief': relief, 'inside': inside,
        'positions': positions, 'is_node': is_node, 'contrast': contrast, 'used': used,
        'box_size': box_size, 'periods': periods, 'image_shape': image.shape,
    }

def xt_correction_to_target(ks, center_rc, target_rc, corr_scale=1e-5):
    """
    X-tilt change that moves the cross from center_rc to target_rc (row, col).
    Moving the fringe pattern by delta_r = target - center changes the phase
    of each fringe set by -2 pi q.delta_r, with q = ks * ronchiPixelSize in
    cycles per pixel. As in ronchi_lib.analyze_ronchigram, the phase error
    (measured - target) = 2 pi q.(target - center) is converted to an xt change
    with ronchi_sem_lib.ronchiCorrMatrix * corr_scale:
    xt x = err_a*corr[0,0] + err_b*corr[0,1], xt y = err_a*corr[1,0] + err_b*corr[1,1].
    The phase error is not wrapped since the cross may be several fringes away.
    Returns (phase_err, xt_delta).
    """
    corr = np.asarray(ronchi_sem_lib.ronchiCorrMatrix, dtype=float) * corr_scale
    q = np.asarray(ks, dtype=float) * ronchi_sem_lib.ronchiPixelSize    # rows: fringe sets a, b
    phase_err = 2 * np.pi * q @ (np.asarray(target_rc, dtype=float) - np.asarray(center_rc, dtype=float))
    return phase_err, corr @ phase_err

def _acquire_and_find_cross(c3, **find_kwargs):
    """Acquire a ronchigram at ronchiC3Offset c3 and find the cross. Returns (image, result)."""
    ronchi_sem_lib.ronchiC3Offset = c3
    image = ronchi_sem_lib.acquire_ronchi_image(ronchi_sem_lib.ronchiStartC3Offset, c3,
                                                sem_acquire_preset='T', pass_label='cross')
    ks, phases = measure_ks_phases(image)
    result = find_cross_center(image, ks, phases, **find_kwargs)
    result.update({'c3_offset': c3, 'ks': ks, 'phases': phases})
    return image, result

def on_cross_correction(c3_offsets=None, target_xy=None, ask=True, show_plot=True, corr_scale=1e-5,
                        min_relief=0.3, **find_kwargs):
    """
    Step through ronchi_sem_lib.ronchiC3Offset in c3_offsets (default -90 to -30
    in steps of 10). At each step, acquire a ronchigram, find the cross center
    and apply the xt correction that moves it to target_xy (x, y in pixels,
    default the image center), so the next step is taken at the corrected
    position. A step is not corrected if its cross center is unreliable: the
    contrast relief (1 - 5th / 95th percentile contrast) is below min_relief,
    which happens at high zoom when the whole image is inside the cross, or
    the fitted center is outside the image. Then acquire a final ronchigram
    at the last C3 offset to show the corrected cross. If ask, the user decides to keep the correction or revert
    to the original xt. xt is also reverted if a step fails.
    ronchi_sem_lib.ronchiC3Offset is restored at the end.
    Returns a dict with the steps, the final result and whether it was kept.
    """
    if c3_offsets is None:
        c3_offsets = np.arange(-90, -29, 10)
    c3_offsets = np.asarray(c3_offsets, dtype=float)
    c3_offset0 = ronchi_sem_lib.ronchiC3Offset
    xt_start = np.array(sem.ReportXLensDeflector(2), dtype=float)
    xt_values = [xt_start]
    steps = []
    kept = False
    try:
        for c3 in c3_offsets:
            image, result = _acquire_and_find_cross(c3, **find_kwargs)
            if target_xy is None:
                target_xy = np.array(result['image_shape'][::-1], dtype=float) / 2
            target_rc = np.asarray(target_xy, dtype=float)[::-1]
            phase_err, xt_delta = xt_correction_to_target(result['ks'], result['center'], target_rc, corr_scale)
            applied = result['relief'] >= min_relief and result['inside']
            xt_new = xt_values[-1] + xt_delta if applied else xt_values[-1]
            if applied:
                sem.SetXLensDeflector(2, xt_new[0], xt_new[1])
            xt_values.append(xt_new)
            result.update({'phase_err': phase_err, 'xt_delta': xt_delta, 'applied': applied})
            steps.append(result)
            offset = result['center'][::-1] - target_xy
            print(f'C3 offset {c3:.1f}: cross offset from target (x, y) ({offset[0]:.1f}, {offset[1]:.1f}) px, '
                  f'relief {result["relief"]:.2f}, phase error (a, b) ({phase_err[0]:.2f}, {phase_err[1]:.2f}) rad, '
                  f'xt change ({xt_delta[0]:.6e}, {xt_delta[1]:.6e}) ' + ('applied' if applied else 'NOT applied'))
            if not applied:
                print(f'WARNING: cross center at C3 offset {c3:.1f} is unreliable '
                      f'(relief {result["relief"]:.2f} < {min_relief} or center outside the image); correction skipped')
            if show_plot:
                display_util.addCrossFit(image, result, subtitle=f'C3 offset {c3:.0f}')
        # final check at the last C3 offset with all corrections applied
        final_image, final = _acquire_and_find_cross(c3_offsets[-1], **find_kwargs)
        final_offset = final['center'][::-1] - target_xy
        xt_final = xt_values[-1]
        print(f'final at C3 offset {c3_offsets[-1]:.1f}: cross offset from target (x, y) '
              f'({final_offset[0]:.1f}, {final_offset[1]:.1f}) px, relief {final["relief"]:.2f}')
        if final['relief'] < min_relief or not final['inside']:
            print(f'WARNING: final cross center is unreliable (relief {final["relief"]:.2f} < {min_relief}); '
                  f'judge from the last corrected steps')
        print(f'xt original ({xt_start[0]:.6e}, {xt_start[1]:.6e}) -> corrected ({xt_final[0]:.6e}, {xt_final[1]:.6e}), '
              f'change ({xt_final[0] - xt_start[0]:.6e}, {xt_final[1] - xt_start[1]:.6e})')
        if show_plot:
            display_util.addCrossFit(final_image, final, subtitle=f'final C3 offset {c3_offsets[-1]:.0f}')
            display_util.addCrossSteps(c3_offsets, [s['center'][::-1] for s in steps], target_xy, xt_values,
                                       final['center'][::-1], [s['applied'] for s in steps])
            display_util.showCrossFit()
            display_util.showCrossSteps()
        kept = not ask or input('Keep the on-cross correction ? (Y/N/y/n)').lower() == 'y'
    finally:
        ronchi_sem_lib.ronchiC3Offset = c3_offset0
        if not kept:
            sem.SetXLensDeflector(2, xt_start[0], xt_start[1])
            print(f'xt reverted to original ({xt_start[0]:.6e}, {xt_start[1]:.6e})')
        else:
            print('on-cross correction kept')
    return {'steps': steps, 'final': final, 'target_xy': target_xy, 'xt_values': np.array(xt_values),
            'kept': kept}

if __name__ == '__main__':
    on_cross_correction()
