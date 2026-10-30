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
    ks (1/um) and phases (rad) of the two fringe sets of an image already
    acquired (ronchi_lib.measure_ks_phases with ronchi_sem_lib settings).
    """
    pixel_size_um = ronchi_sem_lib.ronchiPixelSize if pixel_size_um is None else pixel_size_um
    binning = ronchi_sem_lib.ronchiBinning if binning is None else binning
    peak_radius = ronchi_sem_lib.ronchiPeakRadius if peak_radius is None else peak_radius
    return ronchi_lib.measure_ks_phases(image, pixel_size_um, binning, peak_radius)

def find_cross_center(image, ks, phases, **kwargs):
    """
    Center (row, col) of the laser fringe cross (ronchi_lib.find_cross_center
    with ronchi_sem_lib.ronchiPixelSize and ronchiBinning).
    """
    return ronchi_lib.find_cross_center(image, ks, phases, ronchi_sem_lib.ronchiPixelSize,
                                        ronchi_sem_lib.ronchiBinning, **kwargs)

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
