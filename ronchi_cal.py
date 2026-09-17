#!Python
##############################################################################
# Calibration for Ronchigram  - image analysis with SerialEM calls)
#
##############################################################################
import os
import numpy as np
import sys
import math
import platform
if platform.system() == 'Windows':
    is_simu = False
    sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
else:
    is_simu = True
    print('testing on Mac/Linux with simulator')
import ronchi_lib
import ronchi_sem_lib
import cal_util
import display_util

import serialem as sem

count = 0

def log(text, color=0, style=0):
    if text.startswith("DEBUG:") and not debug:
        return
    if text.startswith("NOTE:"):
        color = 4
    elif text.startswith("WARNING:"):
        color = 5
    elif text.startswith("ERROR:"):
        color = 2
        style = 1
    elif text.startswith("DEBUG:"):
        color = 1
        if breakpoints:
            breakpoint()
    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)

def measure_ronchigram_ks_phases_ronchi_lib(pixel_size_um, binning,
                       peak_radius=100, corr_scale=1e-5, pass_label=''):
    """Set C3Offset, acquire ronchigram and analyze with T preset.
    Use the analysis in ronchi_lib to measure.  This is faster.
    """
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    pass_label = 'cal'
    image = ronchi_sem_lib.acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=pass_label)
    display_util.image_buffer.append(image)
    binned = ronchi_lib._ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = ronchi_lib._ronchi_find_fourier_centered(binned)
    if is_simu:
        import mrcfile
        global count
        for my_tuple in (('ronchi_cal',image),('ronchi_cal_pow', np.abs(image_fft)**2)):
            my_name, my_arr = my_tuple
            with mrcfile.new(f"{my_name}{count:02d}.mrc", overwrite=True) as mrc:
                if np.issubdtype(my_arr.dtype, np.integer) or np.issubdtype(my_arr.dtype, np.float64):
                    mrc.set_data(my_arr.astype(np.float32))
                else:
                    mrc.set_data(my_arr)
        count += 1
    return ronchi_lib._ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius, binning=1,
                                       fourier_size=image_fft.shape[0])

def measure_ronchigram_ks_phases_real_space(pixel_size_um, binning,
                       peak_radius=100, corr_scale=1e-5, pass_label=''):
    """
    Set C3Offset, acquire ronchigram and analyze with T preset
    and then calculate ks and phases using real space fitting.
    This is slower but more accurate.
    """
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    pass_label = 'cal'
    image = ronchi_sem_lib.acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=pass_label)
    display_util.image_buffer.append(image)
    from ronchi_real_space_fit import lppfit
    results = lppfit.run_2d_fringe_fit(image)
    print('xxx leginon real space results')
    ks = np.array([results[1]['wave_freq'],results[2]['wave_freq']])/pixel_size_um
    phases = np.array([results[1]['wave_phase'],results[2]['wave_phase']])*math.pi/180.0
    return ks, phases

def measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius=100, corr_scale=1e-5, pass_label=''):
    return measure_ronchigram_ks_phases_real_space(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)

    #return measure_ronchigram_ks_phases_ronchi_lib(pixel_size_um, binning,
    #                   peak_radius, corr_scale, pass_label)


def calibrate_ronchigram_phase_correction_matrix(pixel_size_um, binning,
                       measure_scope_shift, peak_radius=100, corr_scale=1e-5,
                       c3_correction_factor=20 / 6.85):
    """FFT peak phases -> laser deflector (dict with correction_x/y)."""
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    pass_label = 'cal'
    ks0, phase0 = measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)
    target_phase_a = phase0[0]
    target_phase_b = phase0[1]
    scope_changes = corr_scale * measure_scope_shift * np.array([[1.,0.],[0.,1.]])
    phase_shifts = np.array([[1,0],[0,1]],dtype=float)
    for axis in (0,1):
        start_xlens = np.array((ronchi_sem_lib.ronchiStartXLensX,ronchi_sem_lib.ronchiStartXLensY))
        new_xlens = start_xlens + scope_changes[axis]
        sem.SetXLensDeflector(2, new_xlens[0], new_xlens[1])
        ks, phases = measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)
        sem.SetXLensDeflector(2, start_xlens[0], start_xlens[1])
        # use the definition in analyze_ronchigram
        phase_err_a = np.mod(phases[0] - target_phase_a + np.pi, 2 * np.pi) - np.pi
        phase_err_b = np.mod(phases[1] - target_phase_b + np.pi, 2 * np.pi) - np.pi
        phase_shifts[axis] = np.array([phase_err_a,phase_err_b])
        print(f'ks, phases for axis {axis:d}: {ks}, {phases}')
    print(f'phase_shifts {phase_shifts}')
    if display_util.image_buffer:
        display_util.showImages()
    corr = cal_util.solveTransform(scope_changes, phase_shifts)
    print('correction by xt matrix', corr)
    #global ronchi_sem_lib.ronchiCorrMatrix
    if input('Is this a good matrix ? (Y/N/y/n)').lower() == 'y':
        ronchi_sem_lib.ronchiCorrMatrix = corr.tolist()
        cal_dir, session_name = cal_util.getCalibrationsDir()
        os.makedirs(cal_dir, exist_ok=True)
        cal_util.saveCalibration('ronchi_corr_matrix', cal_dir, session_name,corr) 
    return corr

def correct_phase_error(corr, phase_err_a, phase_err_b):
    correction_x = phase_err_a * corr[0, 0] + phase_err_b * corr[0, 1]
    correction_y = phase_err_a * corr[1, 0] + phase_err_b * corr[1, 1]

def _calibrate_ronchigram_start_c3(c3_delta_scale,
                        pixel_size_um, binning,
                        measure_scope_shift, peak_radius, corr_scale,
                        c3_correction_factor):
    trial_offset_baseline0 = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset0 = ronchi_sem_lib.ronchiC3Offset
    changes = np.array([1.0,0,-1.0]) # n changes
    # (m,n) m sets of n observed shifts
    observed_shifts = np.vstack([changes,changes]) * 0
    observed_shifts = observed_shifts.T
    print('ks',observed_shifts)
    c3d0 = trial_offset_baseline0 + ronchi_offset0
    cal_c3d_changes = c3_delta_scale * changes + ronchi_offset0
    for i, my_change in enumerate(cal_c3d_changes):
        my_c3d = my_change
        ronchi_sem_lib.ronchiC3Offset = my_change
        print(f'c3 change {my_change} in on-plane c3 calibration')
        pass_label = f"my_change:.1f"
        ks0, phase0 = measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)
        my_shift = [np.linalg.norm(ks0[0]),np.linalg.norm(ks0[1])]
        observed_shifts[i] = np.array(my_shift)
    try:
        slopes, intercepts, residuals = cal_util.solveLines(cal_c3d_changes+trial_offset_baseline0, observed_shifts)
        print('intercepts',intercepts)
        print('slopes',slopes)
        print('new start_c3', intercepts/slopes)
        print('residuals', residuals)
    except Exception as e:
        log(f'Error: Calibration not updated {e} Bad pixel shift measured {observed_shifts}')
        return trial_offset_baseline0, np.ones(observed_shifts.shape)*100
    return intercepts/slopes, residuals

def calibrate_ronchigram_start_c3(pixel_size_um, binning,
                       measure_scope_shift, peak_radius=100, corr_scale=1e-5,
                       c3_correction_factor=20 / 6.85):
    """FFT peak ks -> C3imagingdistance saved as the new ronchiStartC3Offset"""
    c3_delta_scale = 5
    c3_offset_diff_threshold = 1
    max_trials = 4
    trial_offset_baseline0 = ronchi_sem_lib.ronchiStartC3Offset
    trial = 1
    while True:
        if trial >= max_trials:
            ronchi_sem_lib.ronchiStartC3Offset = trial_offset_baseline0
            sem.SetImageDistanceOffset(ronchi_sem_lib.ronchiStartC3Offset)
            log('Maximal trials reached. Aborted and c3 offset restored to original')
            return
        c3_offset_arr, residuals = _calibrate_ronchigram_start_c3(c3_delta_scale,
                        pixel_size_um, binning,
                        measure_scope_shift, peak_radius, corr_scale,
                        c3_correction_factor)
        new_trial_offset_baseline = c3_offset_arr.mean()
        if abs(c3_offset_arr[1]-c3_offset_arr[0]) < c3_offset_diff_threshold:
            break
        ronchi_sem_lib.ronchiStartC3Offset -= c3_delta_scale
        sem.SetImageDistanceOffset(ronchi_sem_lib.ronchiStartC3Offset)
        trial += 1
    if display_util.image_buffer:
        display_util.showImages()
    if input('Is this a good ronchiStartC3 ? (Y/N/y/n)').lower() == 'y':
        new_trial_offset_baseline = c3_offset_arr.mean()
        ronchi_sem_lib.ronchiStartC3Offset = new_trial_offset_baseline
        sem.SetImageDistanceOffset(new_trial_offset_baseline)
        cal_dir, session_name = cal_util.getCalibrationsDir()
        os.makedirs(cal_dir, exist_ok=True)
        cal_util.saveCalibration('ronchi_start_c3', cal_dir, session_name,new_trial_offset_baseline) 
        return new_trial_offset_baseline

if __name__=='__main__':
    ronchi_sem_lib.checkRonchigramSetup()
    pixel_size_um = sem.ReportCurrentPixelSize('T')
    xt_tilt = 1 #(scaled at 1e-5 rad)
    ronchi_binning = 32
    corr_scale = 1e-5
    if is_simu:
        corr_scale = 3e-3
        ronchi_binning = 1
        pixel_size_um = 1.5

    # ronchiStartC3 calibration
    new_ronchi_start_c3 = calibrate_ronchigram_start_c3(pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print('new C3Offset',ronchi_sem_lib.ronchiStartC3Offset)

    if False:
        corr_matrix = calibrate_ronchigram_phase_correction_matrix(pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
        print(ronchi_sem_lib.ronchiCorrMatrix)
