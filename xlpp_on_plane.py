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

def refine_c3_by_contrast(step=0.3, max_change=5.0, box_size=512, settle_delay_s=None,
                          apply=True, show_plot=True):
    """
    Refine ImageDistanceOffset (C3) to on-plane by image contrast when xlpp is
    already on-node. Normalized std (std / mean) is measured in a box at the
    image center and in four boxes touching the image corners. Off-plane, the
    corners have lower std/mean than the center. ImageDistanceOffset is scanned
    from start - max_change to start + max_change in increments of step
    (always increasing, to keep lens hysteresis the same), acquiring a Trial
    image at each value, and the offset with the smallest difference between
    the center and the mean of the corners is picked (see _pick_best_c3_offset).
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
    finally:
        # approach the final value from below, like the scan
        sem.SetImageDistanceOffset(offsets[0])
        final_offset = best_offset if apply else start_offset
        sem.SetImageDistanceOffset(final_offset)
        print(f'ImageDistanceOffset set to {final_offset:.2f}')
    if show_plot:
        display_util.addC3ContrastScan(offsets, region_stds, REGION_NAMES, best_offset, start_offset)
    return best_offset, offsets, diffs, region_stds

if __name__=='__main__':
    # ronchiStartC3 calibration
    new_ronchi_start_c3 = ronchi_cal.calibrate_ronchigram_start_c3(ronchi_fixed_pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print('new StartC3Offset',new_ronchi_start_c3)
    # refine C3 to on-plane by center vs corner contrast (xlpp must be on-node)
    best_c3, c3_offsets, c3_diffs, region_stds = refine_c3_by_contrast(step=0.3, max_change=5.0, box_size=512,
                       apply=True, show_plot=True)
    print(f'refined ImageDistanceOffset {best_c3:.2f}')
    if display_util.C3Scan_Offsets:
        display_util.showC3ContrastScan()
