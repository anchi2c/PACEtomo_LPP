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

REGION_NAMES = ('center', 'top-left', 'top-right', 'bottom-left', 'bottom-right')

def _region_normalized_stds(image, box_size=512):
    """
    Normalized standard deviation (std / mean) of five square boxes of the
    image: one at the center and one touching each image corner, in the order
    of REGION_NAMES (numpy row 0 is top). box_size is reduced to half of the
    smaller image dimension if the image is too small for it.
    """
    image = np.asarray(image, dtype=float)
    rows, cols = image.shape
    if box_size > min(rows, cols) // 2:
        box_size = min(rows, cols) // 2
        print(f'WARNING: box size reduced to {box_size} to fit image {rows}x{cols}')
    r_center = rows // 2 - box_size // 2
    c_center = cols // 2 - box_size // 2
    boxes = (
        image[r_center:r_center + box_size, c_center:c_center + box_size],
        image[:box_size, :box_size],
        image[:box_size, cols - box_size:],
        image[rows - box_size:, :box_size],
        image[rows - box_size:, cols - box_size:],
    )
    return np.array([box.std() / box.mean() for box in boxes])

def _pick_best_c3_offset(offsets, diffs, start_offset):
    """
    Offset where center - corners std/mean difference is smallest. If the
    difference changes sign, use the linearly interpolated zero crossing
    nearest to start_offset; otherwise the scanned offset with the smallest
    absolute difference.
    """
    crossings = []
    for k in range(len(offsets) - 1):
        d0, d1 = diffs[k], diffs[k + 1]
        if d0 == 0:
            crossings.append(offsets[k])
        elif d0 * d1 < 0:
            crossings.append(offsets[k] + (offsets[k + 1] - offsets[k]) * d0 / (d0 - d1))
    if diffs[-1] == 0:
        crossings.append(offsets[-1])
    if crossings:
        return float(min(crossings, key=lambda o: abs(o - start_offset))), True
    return float(offsets[np.argmin(np.abs(diffs))]), False

def _fit_flat_cubic(offsets, diffs, low, high, number_of_trials=1001):
    """
    Least-squares fit of a cubic whose inflection is a flat spot,
    diff = a * (offset - x0)**3 + c, with x0 searched on a grid from low to
    high. Returns (x0, coeffs, rss) with coeffs expanded for np.polyval.
    """
    best = None
    for x0 in np.linspace(low, high, number_of_trials):
        design = np.column_stack([(offsets - x0) ** 3, np.ones_like(offsets)])
        (a, c), _, _, _ = np.linalg.lstsq(design, diffs, rcond=None)
        rss = float(np.sum((design @ np.array([a, c]) - diffs) ** 2))
        if best is None or rss < best[2]:
            best = (float(x0), (a, c), rss)
    x0, (a, c), rss = best
    coeffs = np.array([a, -3 * a * x0, 3 * a * x0 ** 2, c - a * x0 ** 3])
    return x0, coeffs, rss

def _pick_best_c3_offset_by_flat_fit(offsets, diffs, center_fraction=0.75, max_rss_ratio=2.0):
    """
    Zero-slope point of a fit of the mean difference (center - mean of corners
    std/mean) against C3 offset, using all scan points. Two models are fitted:
      - a cubic with a flat inflection, a * (offset - x0)**3 + c (flat at x0)
      - a 2nd order polynomial (flat at its vertex)
    A model is used only if its zero-slope point is within the center
    center_fraction of the scanned offset range (the fit is least reliable
    near the ends of the data) and its residual sum of squares is at most
    max_rss_ratio times that of an unconstrained 3rd order polynomial fit (so
    that a model which does not describe the data, e.g. a straight line, is
    not used). If both qualify, the one with the smaller residual sum of
    squares is used (both have 3 parameters).
    Returns (best_offset, coeffs, search_range, info) where coeffs are for
    np.polyval of the chosen fit, search_range is (low, high) and info
    describes the chosen model, or why none was used if best_offset is None.
    """
    offsets = np.asarray(offsets, dtype=float)
    diffs = np.asarray(diffs, dtype=float)
    margin = 0.5 * (1 - center_fraction) * (offsets[-1] - offsets[0])
    low, high = float(offsets[0] + margin), float(offsets[-1] - margin)
    search_range = (low, high)
    if len(offsets) < 4:
        return None, None, search_range, 'fewer than 4 scan points'
    free_coeffs = np.polyfit(offsets, diffs, 3)
    free_rss = float(np.sum((np.polyval(free_coeffs, offsets) - diffs) ** 2))
    max_rss = max_rss_ratio * free_rss
    candidates = []
    # cubic with a flat inflection, constrained to have it within the range
    x0, cubic_coeffs, cubic_rss = _fit_flat_cubic(offsets, diffs, low, high)
    # a flat spot at the edge of the allowed range means the best x0 is outside it
    grid_step = (high - low) / 1000
    if low + grid_step < x0 < high - grid_step and cubic_rss <= max_rss:
        candidates.append((cubic_rss, x0, cubic_coeffs, 'flat-inflection cubic'))
    # 2nd order polynomial
    quad_coeffs = np.polyfit(offsets, diffs, 2)
    quad_rss = float(np.sum((np.polyval(quad_coeffs, offsets) - diffs) ** 2))
    if quad_coeffs[0] != 0:
        vertex = float(-quad_coeffs[1] / (2 * quad_coeffs[0]))
        if low <= vertex <= high and quad_rss <= max_rss:
            candidates.append((quad_rss, vertex, quad_coeffs, '2nd order'))
    if not candidates:
        return None, quad_coeffs, search_range, (
            f'neither fit has its zero slope within the center {center_fraction:.0%} '
            f'of the scan, {low:.2f} to {high:.2f}, with rss <= {max_rss_ratio:g} x '
            f'unconstrained 3rd order {free_rss:.3g}; rss flat-inflection cubic '
            f'{cubic_rss:.3g}, 2nd order {quad_rss:.3g}')
    rss, best_offset, coeffs, model = min(candidates, key=lambda cand: cand[0])
    info = (f'{model} fit, rss flat-inflection cubic {cubic_rss:.3g}, 2nd order {quad_rss:.3g}, '
            f'unconstrained 3rd order {free_rss:.3g}')
    return best_offset, coeffs, search_range, info

def _set_c3_from_below(sem, low_offset, offset):
    """Set ImageDistanceOffset approached from below, like the scan."""
    sem.SetImageDistanceOffset(low_offset)
    sem.SetImageDistanceOffset(offset)

def _ask_user_best_c3(sem, candidates, low_offset, settle_delay_s):
    """
    Acquire a Trial image at each candidate offset and ask the user which is
    better. candidates: list of (label, offset). Returns the chosen offset.
    """
    for label, offset in candidates:
        _set_c3_from_below(sem, low_offset, offset)
        if settle_delay_s > 0:
            sem.Delay(settle_delay_s, "s")
        sem.T()
        print(f'{label}: Trial image acquired at C3 offset {offset:.2f}')
    # each T pushes the previous image down one buffer, so the last is in A
    buffers = 'ABCDEFGH'
    n = len(candidates)
    prompt_lines = []
    for k, (label, offset) in enumerate(candidates):
        prompt_lines.append(f'{k + 1}: {label} C3 offset {offset:.2f} (buffer {buffers[n - 1 - k]})')
    print('\n'.join(prompt_lines))
    while True:
        answer = input('Which image is better? ' + ', '.join(
            f'{k + 1}={label}' for k, (label, _) in enumerate(candidates)) + ' ').strip()
        if answer.isdigit() and 1 <= int(answer) <= n:
            label, offset = candidates[int(answer) - 1]
            print(f'user chose {label} C3 offset {offset:.2f}')
            return offset
        print(f'please enter a number from 1 to {n}')

def refine_c3_by_contrast(step=0.3, max_change=5.0, box_size=512, settle_delay_s=None,
                          apply=True, show_plot=True, ask_user=True):
    """
    Refine ImageDistanceOffset (C3) to on-plane by image contrast when xlpp is
    already on-node. It needs specimen contrast, so it is not valid at a
    vacuum position. Normalized std (std / mean) is measured in a box at the
    image center and in four boxes touching the image corners. Off-plane, the
    corners have lower std/mean than the center. ImageDistanceOffset is scanned
    from start - max_change to start + max_change in increments of step
    (always increasing, to keep lens hysteresis the same), acquiring a Trial
    image at each value, and the offset with the smallest difference between
    the center and the mean of the corners is picked (see _pick_best_c3_offset).
    A second pick is the zero-slope point of a fit of that difference against
    C3 offset, either a cubic with a flat inflection or a 2nd order
    polynomial, within the center 3/4 of the scan
    (see _pick_best_c3_offset_by_flat_fit).
    If ask_user and the two picks differ, a Trial image is acquired at each
    and the user is asked which is better; otherwise the first pick is used.
    The result never moves more than max_change from the starting value.
    settle_delay_s: wait after each change, default ronchi_sem_lib.ronchiDelay.
    apply: set the best offset at the end (approached from below, as in the
    scan); otherwise restore the starting value.
    Returns (best_offset, offsets, diffs, region_stds) where region_stds[k] has
    std/mean of each region in the order of REGION_NAMES at offsets[k].
    """
    sem = ronchi_cal.sem
    if settle_delay_s is None:
        settle_delay_s = ronchi_cal.ronchi_sem_lib.ronchiDelay
    start_offset = float(sem.ReportImageDistanceOffset())
    number_of_steps = int(np.floor(max_change / step + 1e-9))
    offsets = start_offset + step * np.arange(-number_of_steps, number_of_steps + 1)
    region_stds = np.zeros((len(offsets), len(REGION_NAMES)))
    best_offset = start_offset
    try:
        for k, offset in enumerate(offsets):
            sem.SetImageDistanceOffset(offset)
            if settle_delay_s > 0:
                sem.Delay(settle_delay_s, "s")
            sem.T()
            region_stds[k] = _region_normalized_stds(sem.bufferImage("A"), box_size)
            diff = region_stds[k, 0] - region_stds[k, 1:].mean()
            print(f'C3 offset {offset:.2f}: center std/mean {region_stds[k, 0]:.4f}, '
                  f'corners mean {region_stds[k, 1:].mean():.4f}, difference {diff:.4f}')
        diffs = region_stds[:, 0] - region_stds[:, 1:].mean(axis=1)
        best_offset, is_crossing = _pick_best_c3_offset(offsets, diffs, start_offset)
        # never move away from the start by more than max_change
        best_offset = float(np.clip(best_offset, start_offset - max_change, start_offset + max_change))
        how = 'zero crossing' if is_crossing else 'smallest |difference|'
        print(f'best C3 offset {best_offset:.2f} ({how}), start {start_offset:.2f}, change {best_offset - start_offset:+.2f}')
        if not is_crossing and best_offset in (offsets[0], offsets[-1]):
            print(f'WARNING: best C3 offset is at the edge of the +/-{max_change} scan range')
        fit_offset, fit_coeffs, fit_range, fit_info = _pick_best_c3_offset_by_flat_fit(offsets, diffs)
        if fit_offset is None:
            print(f'WARNING: flat fit ignored ({fit_info}); using the mean difference pick only')
        else:
            fit_offset = float(np.clip(fit_offset, start_offset - max_change, start_offset + max_change))
            print(f'flat fit C3 offset {fit_offset:.2f} ({fit_info}), '
                  f'change {fit_offset - start_offset:+.2f}')
        if ask_user and fit_offset is not None and abs(fit_offset - best_offset) > 1e-6:
            best_offset = _ask_user_best_c3(sem, [('mean difference', best_offset),
                                                  ('flat fit', fit_offset)],
                                            offsets[0], settle_delay_s)
        elif ask_user and fit_offset is not None:
            print('both methods picked the same C3 offset')
    finally:
        # approach the final value from below, like the scan
        final_offset = best_offset if apply else start_offset
        _set_c3_from_below(sem, offsets[0], final_offset)
        print(f'ImageDistanceOffset set to {final_offset:.2f}')
    if show_plot:
        display_util.addC3ContrastScan(offsets, region_stds, REGION_NAMES, best_offset, start_offset,
                                       fit_coeffs, fit_offset, fit_range)
    return best_offset, offsets, diffs, region_stds

if __name__=='__main__':
    # ronchiStartC3 calibration
    new_ronchi_start_c3 = ronchi_cal.calibrate_ronchigram_start_c3(ronchi_fixed_pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print('new StartC3Offset',new_ronchi_start_c3)
    # refine C3 to on-plane by center vs corner contrast (xlpp must be on-node).
    # It needs specimen contrast, so it is not valid at a vacuum position.
    if input('Refine on-plane condition using image contrast ? '
             'Only valid on specimen, not at a vacuum position. (Y/N/y/n)').lower() == 'y':
        best_c3, c3_offsets, c3_diffs, region_stds = refine_c3_by_contrast(step=0.3, max_change=5.0, box_size=512,
                           apply=True, show_plot=True)
        print(f'refined ImageDistanceOffset {best_c3:.2f}')
        if display_util.C3Scan_Offsets:
            display_util.showC3ContrastScan()
    else:
        print('on-plane refinement by contrast skipped')
