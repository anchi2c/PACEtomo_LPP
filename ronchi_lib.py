#!Python
##############################################################################
# Ronchigram - image analysis (numpy / FFT only; no SerialEM calls)
#
# Public API: analyze_ronchigram(), measure_ks_phases(), find_cross_center(),
#             fringe_shift()
# Author:        Josh Dickerson, Eric Cooper
# Internal helpers: _ronchi_bin_image, _ronchi_find_fourier_centered, ...
# Extracted from PACEtomo.py by Anchi Cheng
##############################################################################
import numpy as np

def _ronchi_bin_image(image, binning=32):
    bins = [image[i:(image.shape[0] // binning * binning):binning, j:(image.shape[1] // binning * binning):binning]
            for i in range(binning) for j in range(binning)]
    return np.sum(bins, axis=0)


def _ronchi_find_fourier_centered(image, padded_size=4096):
    fourier = np.fft.fftshift(np.fft.fft2(image - np.mean(image), s=(padded_size, padded_size)))
    center = np.array(image.shape) / 2
    grid_y, grid_x = np.indices((padded_size, padded_size), dtype=float)
    shifted_y = grid_y - padded_size / 2
    shifted_x = grid_x - padded_size / 2
    correction_phase_x = np.exp(2j * np.pi * shifted_x / padded_size * center[1])
    correction_phase_y = np.exp(2j * np.pi * shifted_y / padded_size * center[0])
    return fourier * correction_phase_x * correction_phase_y


def _ronchi_report_angles(ks, start_angle=-135):
    return np.mod(np.arctan2(-ks[:, 0], ks[:, 1]) * 180 / np.pi - start_angle, 360) + start_angle


def _ronchi_find_peaks(fourier, radius=100, npeaks=4):
    fourier_abs = np.abs(fourier).copy()
    peak_locations = []
    phases = []
    size = np.shape(fourier_abs)[0]
    peak_coords = [size // 2, size // 2]
    fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0
    for _ in range(npeaks):
        max_idx = np.argmax(fourier_abs)
        peak_coords = np.unravel_index(max_idx, fourier_abs.shape)
        peak_locations.append(peak_coords)
        fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                    (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0
        phases.append(np.angle(fourier[peak_coords]))
    return np.array(peak_locations) - size / 2, np.array(phases)

def analyze_ronchigram_stig(ks):
    """
    Analyze ks in magnitude difference ratio (1= equal) and angle between
    the two vectors.
    """
    norm_0 = np.linalg.norm(ks[0])
    norm_1 = np.linalg.norm(ks[1])
    print(f'magnitudes: {norm_0:.3e}, {norm_1:.3e}')
    norm_diff_ratio = (norm_0-norm_1)/((norm_0+norm_1)/2)
    print(f'magnitude difference: {abs(norm_diff_ratio)*100:.2f} %')
    cos_theta = np.dot(ks[0], ks[1]) / (norm_0 * norm_1)
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_degrees = np.degrees(angle_rad)
    print(f'angle in degrees between the two axis {np.degrees(angle_rad):.2}')
    return norm_diff_ratio, angle_degrees

def _ronchi_find_ks_phases(corrected_fourier, pixel_size_um, npeaks=2, radius=100, binning=1, fourier_size=None):
    """
    image_pixe_size_um: refers to the pixel size of the image acquired by acquire_ronchigram.
    binning: the binning used to shrink the image before calculate corrected_fourier.

    Notes: ks returned is only usable for correction if image_pixel_size_um used in
    magic_phase reference image is at the same binning and image camera_binning and scope
    parameters saved in ronch_ref_ks.jsonl
    """
    if fourier_size is None:
        fourier_size = corrected_fourier.shape[0]
    peaks, phases = _ronchi_find_peaks(corrected_fourier, radius=radius, npeaks=npeaks * 2)
    ordering = np.argsort(_ronchi_report_angles(peaks, start_angle=-135))
    peaks = peaks[ordering][:npeaks]
    phases = phases[ordering][:npeaks]
    ks = peaks / fourier_size * 1 / (pixel_size_um * binning)
    return ks, phases

def analyze_ronchigram(image, pixel_size_um, binning, target_phase_a, target_phase_b,
                       correct_ks, peak_radius=100, corr_matrix=None, corr_scale=1e-5,
                       c3_correction_factor=20 / 6.85):
    """FFT peak phases -> laser deflector and C3 corrections (dict with correction_x/y, phases, ks, ...)."""
    if corr_matrix is None:
        corr_matrix_example = [[0.212, 1.28], [1.22, -0.243]]
        raise ValueError(f'Please define phase to x_tilt correction_matrix at corr_scale of {corr_scale:3f}. For example {corr_matrix_example}')
    correct_ks = np.asarray(correct_ks, dtype=float)
    binned = _ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = _ronchi_find_fourier_centered(binned)
    ks, phases = _ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius, binning=1,
                                       fourier_size=image_fft.shape[0])
    ks_error = ks - correct_ks
    ks_total_err = float(np.linalg.norm(ks_error))
    ks_avg_err = (ks_error[0, 0] + ks_error[1, 1]) / 2.0
    c3_correction = ks_avg_err * c3_correction_factor
    phase_err_a = np.mod(phases[0] - target_phase_a + np.pi, 2 * np.pi) - np.pi
    phase_err_b = np.mod(phases[1] - target_phase_b + np.pi, 2 * np.pi) - np.pi
    corr = np.asarray(corr_matrix, dtype=float) * corr_scale
    correction_x = phase_err_a * corr[0, 0] + phase_err_b * corr[0, 1]
    correction_y = phase_err_a * corr[1, 0] + phase_err_b * corr[1, 1]
    return {
        "ks": ks, "phases": phases, "ks_error": ks_error, "ks_total_err": ks_total_err,
        "ks_avg_err": ks_avg_err,
        "c3_correction": c3_correction,
        "phase_err_a": phase_err_a, "phase_err_b": phase_err_b,
        "correction_x": correction_x, "correction_y": correction_y,
    }

##############################################################################
# Laser fringe cross center (from xlpp_on_cross.py)
##############################################################################
def measure_ks_phases(image, pixel_size_um, binning, peak_radius=100):
    """
    ks (1/um) and phases (rad) of the two fringe sets of an image, the same
    analysis as ronchi_cal.measure_ronchigram_ks_phases.
    pixel_size_um: unbinned pixel size; binning: binning before the FFT.
    Phases are referenced to the center of the binned image.
    """
    binned = _ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = _ronchi_find_fourier_centered(binned)
    return _ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius,
                                  binning=1, fourier_size=image_fft.shape[0])

def _ronchi_phase_reference(image_shape, binning):
    """Unbinned (row, col) position of the binned image center used as phase reference."""
    binned_shape = np.array(image_shape) // binning
    # binned pixel m covers unbinned m*binning .. m*binning + binning - 1
    return binned_shape / 2.0 * binning + (binning - 1) / 2.0

def _ronchi_fringe_lattice(image_shape, ks, phases, pixel_size_um, binning):
    """
    Positions (row, col) in unbinned pixels where both fringe phases are a
    multiple of pi: nodes where both are 0 (mod 2 pi) and antinodes where both
    are pi (mod 2 pi), for the fringe model cos(2 pi k.r + phase).
    Returns (positions, lattice_index (a, b), is_node) for points in the image.
    """
    q = np.asarray(ks, dtype=float) * pixel_size_um   # cycles per unbinned pixel, rows are the two k
    center = _ronchi_phase_reference(image_shape, binning)
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

def _ronchi_lattice_contrast(image, positions, indices, is_node, box_size, beam_box_size, beam_fraction=0.25):
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

def _ronchi_fit_quadratic_2d(x, y, values):
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

def find_cross_center(image, ks, phases, pixel_size_um, binning, box_fraction=0.25, fit_fraction=0.3,
                      beam_fraction=0.25, top_fraction=None, reference_points=None):
    """
    Center (row, col) of the laser fringe cross in unbinned pixels, from the
    local fringe contrast at the nodes and antinodes of the two fringe sets
    and a 2D quadratic fit. ks, phases are from measure_ks_phases with the
    same pixel_size_um and binning.
    box_fraction: box size as a fraction of the shorter fringe period.
    fit_fraction: only points with contrast above fit_fraction of the 95th
    percentile contrast are fitted, so the quadratic describes the cross and
    not the flat area outside it.
    top_fraction: if given, fit only this fraction of the points with the
    highest contrast instead of using fit_fraction. The contrast map is a
    "+": a strong core where both fringe sets overlap and weaker arms of one
    fringe set. When the cross is off center, the arms and weak fringes pull
    the fit; fitting only the strongest points keeps it on the core. Do not
    use it when the core fills the image (high zoom), where the strongest
    points are noise peaks.
    reference_points: if given with top_fraction, top_fraction applies to an
    image with this many valid nodes/antinodes and is scaled by
    reference_points / number of valid points (up to 1). When the core fills
    more of the image the fringes are larger (smaller ks) and there are fewer
    points, so a larger fraction of them is fitted.
    Returns a dict with the center, fit and lattice data.
    """
    image = np.asarray(image, dtype=float)
    q = np.asarray(ks, dtype=float) * pixel_size_um
    periods = 1.0 / np.linalg.norm(q, axis=1)
    box_size = max(3, int(box_fraction * periods.min()))
    positions, indices, is_node = _ronchi_fringe_lattice(image.shape, ks, phases, pixel_size_um, binning)
    box_means, contrast = _ronchi_lattice_contrast(image, positions, indices, is_node, box_size,
                                           int(periods.max()), beam_fraction)
    valid = ~np.isnan(contrast)
    if valid.sum() < 6:
        raise ValueError(f'only {valid.sum()} valid lattice points, can not fit the cross')
    number_of_points = int(valid.sum())
    fraction_used = None
    if top_fraction is not None:
        fraction_used = top_fraction
        if reference_points is not None:
            fraction_used = min(1.0, top_fraction * reference_points / number_of_points)
        # fit only the strongest points, those of the cross core
        threshold = np.percentile(contrast[valid], 100 * (1 - fraction_used))
    else:
        threshold = fit_fraction * np.percentile(contrast[valid], 95)
    used = valid & (contrast >= threshold)
    if used.sum() < 6:
        used = valid
    coeffs, row_max, col_max, is_maximum = _ronchi_fit_quadratic_2d(positions[used, 0], positions[used, 1], contrast[used])
    if not is_maximum:
        print('WARNING: quadratic fit has no maximum, using the highest contrast point')
        row_max, col_max = positions[valid][np.argmax(contrast[valid])]
    # how much the contrast varies across the image; near 0 when the whole
    # image is inside the cross and its center can not be located
    p5, p95 = np.percentile(contrast[valid], [5, 95])
    relief = float(1 - p5 / p95) if p95 > 0 else 0.0
    inside = bool(0 <= row_max < image.shape[0] and 0 <= col_max < image.shape[1])
    # distance of the center from the nearest image edge in fringe periods;
    # the core is partly outside the image when this is small
    edge_distance = min(row_max, col_max, image.shape[0] - 1 - row_max, image.shape[1] - 1 - col_max)
    edge_distance_periods = float(edge_distance / periods.max())
    return {
        'center': np.array([row_max, col_max]), 'is_maximum': is_maximum, 'coeffs': coeffs,
        'relief': relief, 'inside': inside, 'edge_distance_periods': edge_distance_periods,
        'number_of_points': number_of_points, 'fraction_used': fraction_used,
        'positions': positions, 'is_node': is_node, 'contrast': contrast, 'used': used,
        'box_size': box_size, 'periods': periods, 'image_shape': image.shape,
    }

def fringe_shift(ks0, phases0, center0, ks1, phases1, center1, pixel_size_um, max_ks_change=0.05):
    """
    Shift (row, col) in unbinned pixels of the fringe pattern from image 0 to
    image 1 of the same shape. A shift d changes the phase of each fringe set
    by -2 pi q.d (q = ks * pixel_size_um, cycles per pixel), which gives d up
    to whole fringes. The number of whole fringes is the one closest to the
    cross center shift center1 - center0 (from find_cross_center), so the
    cross center error must be less than half a fringe period.
    If ks change by more than max_ks_change (fraction), the pattern changed
    by more than a shift and the cross center shift is returned instead.
    Returns (shift, info) where info has 'center_shift', 'fringes' (whole
    fringes added to each phase), 'ambiguity' (distance of the center-predicted
    fringe count from the nearest whole number, 0.5 is fully ambiguous),
    'ks_change' and 'phase_refined'.
    """
    ks0 = np.asarray(ks0, dtype=float)
    ks1 = np.asarray(ks1, dtype=float)
    center_shift = np.asarray(center1, dtype=float) - np.asarray(center0, dtype=float)
    ks_change = float(np.max(np.linalg.norm(ks1 - ks0, axis=1) / np.linalg.norm(ks0, axis=1)))
    info = {'center_shift': center_shift, 'ks_change': ks_change, 'fringes': None, 'ambiguity': None,
            'phase_refined': False}
    if ks_change > max_ks_change:
        return center_shift, info
    q = (ks0 + ks1) / 2.0 * pixel_size_um
    dphi = np.asarray(phases1, dtype=float) - np.asarray(phases0, dtype=float)
    # whole fringes that make the phase change agree with the cross center shift
    fringe_count = (-2 * np.pi * q @ center_shift - dphi) / (2 * np.pi)
    fringes = np.round(fringe_count)
    shift = np.linalg.solve(-2 * np.pi * q, dphi + 2 * np.pi * fringes)
    info.update({'fringes': fringes, 'ambiguity': float(np.max(np.abs(fringe_count - fringes))),
                 'phase_refined': True})
    return shift, info
