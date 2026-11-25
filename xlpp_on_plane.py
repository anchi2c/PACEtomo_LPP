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

def power_inflection(offsets, x0, power, a, c):
    """
    a * sign(offset - x0) * |offset - x0|**power + c. For power > 1 the slope
    is zero at x0, and the flat part around x0 is wider for a larger power.
    """
    u = np.asarray(offsets, dtype=float) - x0
    return a * np.sign(u) * np.abs(u) ** power + c

def _fit_power_inflection(offsets, diffs, low, high, powers, number_of_x0=1001):
    """
    Least-squares fit of power_inflection to diffs, with x0 searched on a grid
    from low to high and power over powers; a and c are solved linearly for
    each (x0, power). Returns (params, rss) with params a dict of
    x0, power, a, c.
    """
    x0s = np.linspace(low, high, number_of_x0)
    y = diffs - diffs.mean()
    best = None
    for power in powers:
        u = offsets[np.newaxis, :] - x0s[:, np.newaxis]
        basis = np.sign(u) * np.abs(u) ** power
        basis_centered = basis - basis.mean(axis=1, keepdims=True)
        var = np.sum(basis_centered ** 2, axis=1)
        cov = basis_centered @ y
        valid = var > 0
        rss = np.full(len(x0s), np.inf)
        rss[valid] = np.sum(y ** 2) - cov[valid] ** 2 / var[valid]
        k = int(np.argmin(rss))
        if best is None or rss[k] < best[1]:
            a = cov[k] / var[k]
            c = diffs.mean() - a * basis[k].mean()
            best = ({'x0': float(x0s[k]), 'power': float(power), 'a': float(a), 'c': float(c)},
                    float(max(rss[k], 0.0)))
    return best

def _pick_best_c3_offset_by_flat_fit(offsets, diffs, center_fraction=0.75, max_rss_ratio=2.0,
                                     max_line_rss_ratio=0.5, powers=np.arange(1.2, 8.01, 0.1)):
    """
    Flat point of a fit of the mean difference (center - mean of corners
    std/mean) against C3 offset with a general power inflection function,
    a * sign(offset - x0) * |offset - x0|**power + c (see power_inflection),
    using all scan points. The slope is zero at x0, and the power sets how
    wide the flat part is. x0 is searched only within the center
    center_fraction of the scanned offset range (the fit is least reliable
    near the ends of the data), and power over powers.
    The fit is used only if x0 is not at the edge of that range, and its
    residual sum of squares is at most max_rss_ratio times that of an
    unconstrained 3rd order polynomial fit (so that the fit is not used for
    data it does not describe) and at most max_line_rss_ratio times that of a
    straight line fit (so that data without a flat part, where a power near
    1 fits as well as a line with an arbitrary x0, is not used).
    Returns (best_offset, params, search_range, info) where params is a dict
    of x0, power, a, c (None if fewer than 4 points), search_range is
    (low, high) and info describes the fit, or why it was not used if
    best_offset is None.
    """
    offsets = np.asarray(offsets, dtype=float)
    diffs = np.asarray(diffs, dtype=float)
    margin = 0.5 * (1 - center_fraction) * (offsets[-1] - offsets[0])
    low, high = float(offsets[0] + margin), float(offsets[-1] - margin)
    search_range = (low, high)
    if len(offsets) < 4:
        return None, None, search_range, 'fewer than 4 scan points'
    params, rss = _fit_power_inflection(offsets, diffs, low, high, powers)
    free_coeffs = np.polyfit(offsets, diffs, 3)
    free_rss = float(np.sum((np.polyval(free_coeffs, offsets) - diffs) ** 2))
    line_coeffs = np.polyfit(offsets, diffs, 1)
    line_rss = float(np.sum((np.polyval(line_coeffs, offsets) - diffs) ** 2))
    info = (f'power {params["power"]:.1f}, rss {rss:.3g}, '
            f'unconstrained 3rd order {free_rss:.3g}, line {line_rss:.3g}')
    # a best x0 at the edge of the grid means the flat point is outside the range
    grid_step = (high - low) / 1000
    if not low + grid_step < params['x0'] < high - grid_step:
        return None, params, search_range, (
            f'flat point at the edge of the center {center_fraction:.0%} of the scan, '
            f'{low:.2f} to {high:.2f}; {info}')
    if rss > max_rss_ratio * free_rss:
        return None, params, search_range, (
            f'rss > {max_rss_ratio:g} x unconstrained 3rd order; {info}')
    if rss > max_line_rss_ratio * line_rss:
        return None, params, search_range, (
            f'not clearly better than a straight line, rss > {max_line_rss_ratio:g} x line; {info}')
    return params['x0'], params, search_range, info

def _set_c3_from_below(sem, low_offset, offset):
    """Set ImageDistanceOffset approached from below, like the scan."""
    sem.SetImageDistanceOffset(low_offset)
    sem.SetImageDistanceOffset(offset)

def refine_c3_by_contrast(step=0.3, max_change=5.0, box_size=512, settle_delay_s=None,
                          apply=True, show_plot=True, fit_offset_correction=-1.0):
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
    A second pick is the flat point of a power inflection function fit of
    that difference against C3 offset, within the center 3/4 of the scan
    (see _pick_best_c3_offset_by_flat_fit).
    The flat-fit pick plus fit_offset_correction is used when the fit is valid
    (the fit pick is usually at the high end, so the default correction is
    -1.0); otherwise the first pick is used.
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
        print(f'mean difference C3 offset {best_offset:.2f} ({how}), start {start_offset:.2f}, change {best_offset - start_offset:+.2f}')
        if not is_crossing and best_offset in (offsets[0], offsets[-1]):
            print(f'WARNING: mean difference C3 offset is at the edge of the +/-{max_change} scan range')
        fit_offset, fit_params, fit_range, fit_info = _pick_best_c3_offset_by_flat_fit(offsets, diffs)
        if fit_offset is None:
            print(f'WARNING: flat fit ignored ({fit_info}); using the mean difference pick only')
        else:
            fit_offset = float(np.clip(fit_offset, start_offset - max_change, start_offset + max_change))
            print(f'flat fit C3 offset {fit_offset:.2f} ({fit_info}), '
                  f'change {fit_offset - start_offset:+.2f}')
            # the flat-fit pick is usually at the high end
            best_offset = float(np.clip(fit_offset + fit_offset_correction,
                                        start_offset - max_change, start_offset + max_change))
            print(f'using flat fit C3 offset {fit_offset:.2f} {fit_offset_correction:+.2f} = {best_offset:.2f}, '
                  f'change {best_offset - start_offset:+.2f}')
    finally:
        # approach the final value from below, like the scan
        final_offset = best_offset if apply else start_offset
        _set_c3_from_below(sem, offsets[0], final_offset)
        print(f'ImageDistanceOffset set to {final_offset:.2f}')
    if show_plot:
        fit_curve = fit_label = None
        if fit_params is not None:
            fine = np.linspace(offsets[0], offsets[-1], 200)
            fit_curve = (fine, power_inflection(fine, fit_params['x0'], fit_params['power'],
                                                fit_params['a'], fit_params['c']))
            fit_label = f'power inflection fit (power {fit_params["power"]:.1f})'
        display_util.addC3ContrastScan(offsets, region_stds, REGION_NAMES, best_offset, start_offset,
                                       fit_curve, fit_offset, fit_range, fit_label)
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
