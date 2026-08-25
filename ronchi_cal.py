#!Python
##############################################################################
# Calibration for Ronchigram  - image analysis with SerialEM calls)
#
##############################################################################
import os
import numpy as np
import ronchi_lib
import ronchi_sem_lib
import cal_util
import serialem as sem

is_simu = False
count = 0

def measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius=100, corr_scale=1e-5, pass_label=''):
    """Set C3Offset, acquire ronchigram and analyze with T preset"""
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    pass_label = 'cal'
    image = ronchi_sem_lib.acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=pass_label)
    binned = ronchi_lib._ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = ronchi_lib._ronchi_find_fourier_centered(binned)
    if is_simu:
        global count
        from pyami import mrc
        mrc.write(image,f'image{count:d}.mrc')
        mrc.write(np.abs(image_fft),f'power{count}.mrc')
        count += 1
    return ronchi_lib._ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius, binning=1,
                                       fourier_size=image_fft.shape[0])

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
    corr = cal_util.solveTransform(scope_changes, phase_shifts)
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

def calibrate_ronchigram_start_c3(pixel_size_um, binning,
                       measure_scope_shift, peak_radius=100, corr_scale=1e-5,
                       c3_correction_factor=20 / 6.85):
    """FFT peak ks -> C3imagingdistance saved as the new ronchiStartC3Offset"""
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset0 = ronchi_sem_lib.ronchiC3Offset
    pass_label = '_cal0'
    ks0, phase0 = measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)
    delta_offset = 0.5 * ronchi_sem_lib.ronchiC3Offset
    ronchi_sem_lib.ronchiC3Offset += delta_offset
    pass_label = '_cal1'
    ks1, phase1 = measure_ronchigram_ks_phases(pixel_size_um, binning,
                       peak_radius, corr_scale, pass_label)
    ronchi_sem_lib.ronchiC3Offset = ronchi_offset0
    # fit as: ks = ronchiStartC3Offset + a * ronchiC3Offset
    # need to solve for a, the slope and betterronchiStartC3Offset
    a = (ks1 - ks0)/delta_offset
    new_offset_baseline_array = ks0 - a*ronchi_offset0
    lengths = np.linalg.norm(new_offset_baseline_array, axis=1)
    new_trial_offset_baseline = np.mean(lengths)
    print(new_trial_offset_baseline)
    if input('Is this a good ronchiStartC3 ? (Y/N/y/n)').lower() == 'y':
        ronchi_sem_lib.ronchiStartC3 = new_trial_offset_baseline
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
    
    corr_matrix = calibrate_ronchigram_phase_correction_matrix(pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print(ronchi_sem_lib.ronchiCorrMatrix)
 
    # ronchiStartC3 calibration
    new_ronchi_start_c3 = calibrate_ronchigram_start_c3(pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print(ronchi_sem_lib.ronchiStartC3)
