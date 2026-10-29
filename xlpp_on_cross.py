#!Python
##############################################################################
# On-cross measurement and correction for xlpp.
#
# The two laser fringe sets cross in the ronchigram. The fringe contrast is
# highest at the center of the cross. Changing the ronchigram C3 offset zooms
# the ronchigram about the beam axis by about 1/|C3 offset|, so a cross that
# is not on the axis moves across the scan. On-cross, it does not move.
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
    return {
        'center': np.array([row_max, col_max]), 'is_maximum': is_maximum, 'coeffs': coeffs,
        'positions': positions, 'is_node': is_node, 'contrast': contrast, 'used': used,
        'box_size': box_size, 'periods': periods, 'image_shape': image.shape,
    }

def measure_on_cross(c3_offsets=None, show_plot=True, **find_kwargs):
    """
    Acquire a ronchigram at each ronchi_sem_lib.ronchiC3Offset in c3_offsets
    (default -90 to -30 in steps of 10), measure ks and phases, and find the
    center of the laser fringe cross. Then fit center = axis + drift / |C3 offset|
    for x and y. axis is the zoom center (beam axis) and drift is zero on-cross.
    The fit is weighted by |C3 offset| because the center error grows with the
    zoom (about 1/|C3 offset|), and at high zoom more of the cross is outside
    the image.
    ronchi_sem_lib.ronchiC3Offset is restored at the end.
    Returns a dict with per-scan results and the fit.
    """
    if c3_offsets is None:
        c3_offsets = np.arange(-90, -29, 10)
    c3_offsets = np.asarray(c3_offsets, dtype=float)
    baseline = ronchi_sem_lib.ronchiStartC3Offset
    c3_offset0 = ronchi_sem_lib.ronchiC3Offset
    scans = []
    try:
        for c3 in c3_offsets:
            ronchi_sem_lib.ronchiC3Offset = c3
            image = ronchi_sem_lib.acquire_ronchi_image(baseline, c3, sem_acquire_preset='T', pass_label='cross')
            ks, phases = measure_ks_phases(image)
            result = find_cross_center(image, ks, phases, **find_kwargs)
            result.update({'c3_offset': c3, 'ks': ks, 'phases': phases})
            scans.append(result)
            row, col = result['center']
            print(f'C3 offset {c3:.1f}: fringe periods {np.array2string(result["periods"], precision=1)} px, '
                  f'{int(result["used"].sum())} points fitted, cross center (x, y) = ({col:.1f}, {row:.1f})')
            if show_plot:
                display_util.addCrossFit(image, result, subtitle=f'C3 offset {c3:.0f}')
    finally:
        ronchi_sem_lib.ronchiC3Offset = c3_offset0
    inv_c3 = 1.0 / np.abs(c3_offsets)
    centers_xy = np.array([s['center'][::-1] for s in scans])    # sem convention (x, y)
    # center = axis + drift * (1/|c3|), fitted for x and y, weighted by |c3|
    weights = np.abs(c3_offsets) / np.abs(c3_offsets).max()
    design = np.column_stack([np.ones_like(inv_c3), inv_c3]) * weights[:, None]
    (axis_xy, drift_xy), *_ = np.linalg.lstsq(design, centers_xy * weights[:, None], rcond=None)
    print(f'beam axis (zoom center) (x, y) = ({axis_xy[0]:.1f}, {axis_xy[1]:.1f}) px')
    print(f'cross drift (x, y) = ({drift_xy[0]:.1f}, {drift_xy[1]:.1f}) px*um; '
          f'cross offset from axis at C3 offset {c3_offsets[-1]:.0f}: '
          f'({drift_xy[0] * inv_c3[-1]:.1f}, {drift_xy[1] * inv_c3[-1]:.1f}) px')
    if show_plot:
        display_util.addCrossCenters(c3_offsets, centers_xy, axis_xy, drift_xy)
    return {'scans': scans, 'c3_offsets': c3_offsets, 'centers_xy': centers_xy,
            'axis_xy': axis_xy, 'drift_xy': drift_xy}

def correct_on_cross(measurement, ask=True, apply=True, show_plot=True, corr_scale=1e-5):
    """
    X-tilt change that moves the cross onto the beam axis, extrapolated to
    C3 offset 0 (on-plane).
    For each scan, moving the fringe pattern by delta_r = axis - cross center
    changes the phase of each fringe set by -2 pi q.delta_r, with q = ks *
    ronchiPixelSize in cycles per pixel. As in ronchi_lib.analyze_ronchigram,
    the phase error (measured - target) = 2 pi q.(axis - center) is converted
    to an xt change with ronchi_sem_lib.ronchiCorrMatrix * corr_scale:
    xt change x = err_a*corr[0,0] + err_b*corr[0,1], y = err_a*corr[1,0] + err_b*corr[1,1].
    The phase error is not wrapped since the cross may be several fringes away.
    The xt changes are fitted linearly against C3 offset and the value at
    C3 offset 0 is applied.
    ask: confirm before applying. Returns the xt change.
    """
    corr = np.asarray(ronchi_sem_lib.ronchiCorrMatrix, dtype=float) * corr_scale
    c3_offsets = np.asarray(measurement['c3_offsets'], dtype=float)
    axis_rc = np.asarray(measurement['axis_xy'], dtype=float)[::-1]    # (row, col)
    scan_phase_errs = []
    scan_xt_deltas = []
    for scan in measurement['scans']:
        q = np.asarray(scan['ks'], dtype=float) * ronchi_sem_lib.ronchiPixelSize   # rows: fringe sets a, b
        phase_err = 2 * np.pi * q @ (axis_rc - scan['center'])
        scan_phase_errs.append(phase_err)
        scan_xt_deltas.append(corr @ phase_err)
    scan_phase_errs = np.array(scan_phase_errs)
    scan_xt_deltas = np.array(scan_xt_deltas)
    # linear fit xt change = a + b * c3 offset, for xt x and y; a is the value at C3 offset 0
    design = np.column_stack([np.ones_like(c3_offsets), c3_offsets])
    (xt_delta, xt_slope), *_ = np.linalg.lstsq(design, scan_xt_deltas, rcond=None)
    for c3, e, d in zip(c3_offsets, scan_phase_errs, scan_xt_deltas):
        print(f'C3 offset {c3:.1f}: phase error (a, b) ({e[0]:.2f}, {e[1]:.2f}) rad, xt change ({d[0]:.6e}, {d[1]:.6e})')
    print(f'xt change extrapolated to C3 offset 0: ({xt_delta[0]:.6e}, {xt_delta[1]:.6e}), '
          f'slope ({xt_slope[0]:.3e}, {xt_slope[1]:.3e}) per C3 unit')
    if show_plot:
        display_util.addCrossCorrection(c3_offsets, scan_xt_deltas, xt_delta, xt_slope)
    xt0 = np.array(sem.ReportXLensDeflector(2), dtype=float)
    xt1 = xt0 + xt_delta
    print(f'xt ({xt0[0]:.6e}, {xt0[1]:.6e}) -> ({xt1[0]:.6e}, {xt1[1]:.6e})')
    if not apply:
        return xt_delta
    if ask and input('Apply this on-cross correction ? (Y/N/y/n)').lower() != 'y':
        print('on-cross correction not applied')
        return xt_delta
    sem.SetXLensDeflector(2, xt1[0], xt1[1])
    print('on-cross correction applied')
    return xt_delta

if __name__ == '__main__':
    measurement = measure_on_cross()
    if display_util.Cross_Images:
        display_util.showCrossFit()
    if display_util.CrossCenter_C3:
        display_util.showCrossCenters()
    correct_on_cross(measurement, show_plot=True)
    if display_util.CrossCorrection_C3:
        display_util.showCrossCorrection()
