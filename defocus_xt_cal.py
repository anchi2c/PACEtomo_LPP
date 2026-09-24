#!Python
# ===================================================================
#ScriptName     Defocus XLens Calibration
# Purpose:      Runs Calibration of XLens deflector correction for defocus change.
# Author:       Anchi Cheng
# ===================================================================

############ SETTINGS ############

########## Ronchigram settings ##########
# Requires ronchi_sem_lib.hasXLens = True.
# Settings set in ronchi_sem_lib.py
########## END Ronchigram settings ##########

import sys
import platform
if platform.system() == 'Windows':
    is_simu = False
    sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
    import serialem as sem
else:
    is_simu = True
    print('testing on Mac/Linux with simulator')
    import sem_simulator as sem
import os
import numpy as np
import ronchi_sem_lib
import cal_util
import display_util
import lafis_cal

######### defocus xt correction #########
# Defocus is a scalar while xt is [x,y]. The calibration is a vector
# of xt change per um of defocus change, used as
# xt1 = xt0 + df_delta * df_xt_vector
df_xt_vector = [0.0, 0.0]   #rad/um, zero start makes the first refinement a raw measurement

pixel_xt_matrix = [[1.0866279077999436e-06, -1.7905014130717222e-07], [2.861300907326552e-07, 1.8702918955272723e-06]] #starting guess mrad/pixel

cal_df_scales = [2.5, 5, 10]    #in um. Changes of +/- scale are used. Keep within +/- 10 um
max_trials = 3
converging_deviation_threshold = 15   #in pixels
##########

log = lafis_cal.log

def saveCalibrations():
    cal_dir, session_name = cal_util.getCalibrationsDir()
    os.makedirs(cal_dir, exist_ok=True)
    mag,*_ = sem.ReportMag()
    cal_util.saveCalibration('pixel_xt_matrix_%d' % int(mag), cal_dir, session_name,pixel_xt_matrix)
    cal_util.saveCalibration('df_xt_vector', cal_dir, session_name,df_xt_vector)

def readCalibrations():
    global pixel_xt_matrix, df_xt_vector
    # lafis_cal is_xt_matrix is used to scale xt change in calibrateXtPixelMatrix
    lafis_cal.readCalibrations()
    cal_dir, session_name = cal_util.getCalibrationsDir()
    mag,*_ = sem.ReportMag()
    # read calibrations from file. Only replace hardcoded default if
    # there is saved value
    r = cal_util.readCalibration('pixel_xt_matrix_%d' % int(mag), cal_dir)
    if r:
        pixel_xt_matrix = r
    r = cal_util.readCalibration('df_xt_vector', cal_dir)
    if r:
        df_xt_vector = r

def calc_xt_df(xt0, df_delta):
    # This form works for both array and list [x,y]
    xt1 = [0.0,0.0]
    xt1[0] = xt0[0]+df_delta*df_xt_vector[0]
    xt1[1] = xt0[1]+df_delta*df_xt_vector[1]
    return xt1

def calibrateXtPixelMatrix():
    global pixel_xt_matrix
    lafis_cal.pixel_xt_matrix = pixel_xt_matrix
    lafis_cal.calibrateXtPixelMatrix()
    pixel_xt_matrix = lafis_cal.pixel_xt_matrix

def _measureDefocusXtResidual(cal_df_changes, trial_offset_baseline, ronchi_offset):
    df0 = float(sem.ReportDefocus())
    xt0 = [float(v) for v in sem.ReportXLensDeflector(2)[:2]]
    img0_array = lafis_cal._acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
    residual_shifts = np.zeros((len(cal_df_changes),2), dtype=np.float32)
    applied_xts = []
    for i, my_df in enumerate(cal_df_changes):
        sem.SetDefocus(df0+my_df)
        xt1 = calc_xt_df(xt0, my_df)
        sem.SetXLensDeflector(2, xt1[0], xt1[1])
        img_array = lafis_cal._acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
        display_util.addImage(img_array)
        residual = lafis_cal._find_shift_sem(img0_array, img_array)
        residual_shifts[i] = np.array(residual)
        applied_xts.append([xt1[0]-xt0[0],xt1[1]-xt0[1]])
        sem.SetDefocus(df0)
        sem.SetXLensDeflector(2, xt0[0], xt0[1])
    print('defocus xt residual pixel shift x,y', residual_shifts)
    return np.array(applied_xts), residual_shifts

def update_df_xt_vector(slopes):
    global df_xt_vector
    print('before update df_xt_vector', df_xt_vector)
    df_xt_vector = slopes.tolist()
    print('updated df_xt_vector', df_xt_vector)

def _refineDefocusXtVector(df_scale, trial_offset_baseline, ronchi_offset):
    cal_df_changes = df_scale * np.array([-1.0,-0.5,0.5,1.0])
    print('defocus changes for cal', cal_df_changes)
    applied_xts, pixel_residuals = _measureDefocusXtResidual(cal_df_changes, trial_offset_baseline, ronchi_offset)
    xt_residuals_arr = pixel_residuals @ np.array(pixel_xt_matrix)
    print('applied_xts', applied_xts)
    print('xt_residuals', xt_residuals_arr)
    try:
        # solve x and y of xt against scalar defocus change as two lines
        slopes, intercepts, residuals = cal_util.solveLines(cal_df_changes, applied_xts-xt_residuals_arr)
        print('intercepts', intercepts)
        if np.any(np.abs(intercepts) > 0.1 * np.abs(slopes * df_scale)):
            log(f'WARNING: large intercepts {intercepts} relative to xt change {slopes*df_scale}')
        update_df_xt_vector(slopes)
    except Exception as e:
        log(f'Error: Calibration not updated {e} Bad xt residuals {xt_residuals_arr}')
    return pixel_residuals

def calibrateDefocusXt():
    lafis_cal.checkRonchigramSetup()
    lafis_cal.saveZeroImageShiftDefocusXLens()
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    for scale in cal_df_scales:
        # do a refinement of the existing df_xt_vector
        log(f'calibrating df_xt_vector with defocus change of +/-{scale} um')
        trial = 1
        while True:
            if trial > max_trials:
                raise ValueError('Defocus xt calibration did not converge.')
            log(f'calibrating df_xt_vector trial {trial:d} with defocus change of +/-{scale} um')
            pixel_residual = _refineDefocusXtVector(scale, trial_offset_baseline, ronchi_offset)
            mean_deviation = np.linalg.norm(pixel_residual, axis=1).mean()
            if mean_deviation < converging_deviation_threshold:
                break
            trial += 1
    lafis_cal.resetOptics()

def testDefocusXt():
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    df0 = float(sem.ReportDefocus())
    xt0 = [float(v) for v in sem.ReportXLensDeflector(2)[:2]]
    for my_df in (0, 5, -5, 10, -10):
        sem.SetDefocus(df0+my_df)
        xt1 = calc_xt_df(xt0, my_df)
        sem.SetXLensDeflector(2, xt1[0], xt1[1])
        img1_array = lafis_cal._acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
        display_util.addImage(img1_array)
    lafis_cal.resetOptics()

if __name__=='__main__':
    lafis_cal.checkRonchigramSetup()
    lafis_cal.saveZeroImageShiftDefocusXLens()
    readCalibrations()
    ##### calibrate pixel_xt_matrix
    failed_xt_pixel = False
    try:
        calibrateXtPixelMatrix()
        saveCalibrations()
        if display_util.image_buffer:
            display_util.showImages()
            # reset image_buffer
            display_util.image_buffer = []
    except Exception as e:
        log(f'Failed: {e}')
        log('Aborting....')
        failed_xt_pixel = True

    if not failed_xt_pixel:
        ##### defocus xt
        try:
            calibrateDefocusXt()
            saveCalibrations()
        except Exception as e:
            log(f'Failed: {e}')
            log('Aborting....')
        #testDefocusXt()
    print(f'final pixel_xt_matrix: {pixel_xt_matrix}')
    if not failed_xt_pixel:
        print(f'final df_xt_vector: {df_xt_vector}')
    if display_util.image_buffer:
        display_util.showImages()
