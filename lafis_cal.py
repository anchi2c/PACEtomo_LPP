#!Python
# ===================================================================
#ScriptName     Lafis Calibration
# Purpose:      Runs Calibration for LAFIS.
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
import copy
import time
import struct
import platform
from datetime import datetime, timezone
import json
import glob
import numpy as np
from scipy import optimize, ndimage
import ronchi_sem_lib
import cal_util
import display_util

######### LAFIS: lpp afis correction #########
# calibration matrix applied when beamTiltComp == True on xlpp
# Requires ronchi_sem_lib.hasXLens = True and beamTiltComp = True to be meaningful.
count = 0
# all transform matrix are to be used with [x,y] as in SerialEM convention
is_xt_matrix = [[0.00028763286791300477, -0.00012083706524300338], [0.0009948701965659995, 0.00014791197022821008]]  #simu_testing
is_xt_matrix = [[0.00024180881486977201, -0.00042928950183638887], [0.0010020323076790395, 0.0002672870290525659]]
df_is_matrix = [[0.041381,0.012342], [0.041381,0.012342]]

pixel_xt_matrix = [[1.0866279077999436e-06, -1.7905014130717222e-07], [2.861300907326552e-07, 1.8702918955272723e-06]] #starting guess mrad/pixel

lafisZeroImageShiftDefocus = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisZeroImageShiftXLens = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisIsDone = False            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisXtCorrectionX = 0.0       # set from doLafis as the correction made on XLens
lafisXtCorrectionY = 0.0       # set from doLafis as the correction made on XLens
########## 

# calibration parameters
working_dir = sem.ReportDirectory()
cal_dir, session_name = cal_util.getCalibrationsDir(working_dir)
os.makedirs(cal_dir, exist_ok=True)
optics = ronchi_sem_lib.getOpticsKey(ronchi_binning=None, ronchi_c3_offset=ronchi_sem_lib.ronchiC3Offset)
print(optics)

# ScriptName Script 11 Recall xt0 and other origin values from temp_xt0.json

def getResetOpticsPath():
    if platform.system() == 'Windows':
        # TODO: should use working directory
        working_dir = sem.ReportDirectory()
        filepath = os.path.join(working_dir,'temp_xt0.json')
        
    else:
        filepath = './temp_xt0.json'
    if not os.path.exists(filepath):
        raise ValueError(f'optics not saved in {filepath} to be used for reset')
    return filepath

def resetOptics():
    filepath = getResetOpticsPath()
    sem.Echo('-------- Loading optical values from %s' % os.path.join(os.getcwd(), filepath))

    with open(filepath, "r") as f:
        data = json.load(f)
    sem.SetImageShift(data['image_shift'][0], data['image_shift'][1])
    sem.SetBeamTilt(data['beam_tilt'][0],data['beam_tilt'][1])
    sem.SetXLensDeflector(2, data['x_tilt'][0],data['x_tilt'][1])
    sem.SetObjectiveStigmator(data['obj_stig'][0],data['obj_stig'][1])
    sem.SetDefocus(data['defocus'])
    sem.SetImageDistanceOffset(data['c3_offset'])
    sem.Echo('Value reset')

def breakpoint():
    """Breakpoint for debugging in SerialEM."""
    while not sem.KeyBreak():
        sem.Delay(0.1, "s")
    for i in range(5):
        if sem.KeyBreak("d"):
            dumpVars()
            break
        sem.Delay(0.1, "s")


def log(text, color=0, style=0):
    if text.startswith("DEBUG:") and not debug:
        return
    if text.startswith("NOTE:"):
        color = 4
    elif text.startswith("WARNING:"):
        color = 5
    elif text.startswith("ERROR:"):
        print(text)
        color = 2
        style = 1 
    elif text.startswith("DEBUG:"):
        color = 1
        if breakpoints:
            breakpoint()
    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)

def suggestRonchiC3Offset(optics):
    optics['ronchi_c3_offset']=None
    r, ronchi_c3_offset = cal_util.readCalibration('is_xt_matrix', cal_dir, optics)
    if r:
        return ronchi_c3_offset
    if optics['mag'] == 54000:
        return -125.0 # xt_pixel xt_is 54000 1.35 um spot 6
    elif optics['mag'] in (88000, 110000):
        return -130.0 # xt_pixel xt_is 54000 1.35 um spot 6
    else:
        return -abs(float(input('Enter an C3 offset that gives strong lpp cross')))

def checkRonchigramSetup():
    global optics
    filepath = getResetOpticsPath()
    ronchi_sem_lib.checkRonchigramSetup()
    ronchi_sem_lib.ronchiC3Offset = suggestRonchiC3Offset(optics)
    optics = ronchi_sem_lib.getOpticsKey(ronchi_binning=None, ronchi_c3_offset=ronchi_sem_lib.ronchiC3Offset)
    #ronchi_sem_lib.ronchiC3Offset = -173.0 # xt_pixel xt_is 88000 1.5 um
    #ronchi_sem_lib.ronchiC3Offset = -100.0 # xt_pixel xt_is 54000 1.5 um
    #ronchi_sem_lib.ronchiC3Offset = -90.0 # xt_pixel xt_is 54000 1.35 um spot 6
    #ronchi_sem_lib.ronchiC3Offset = -130.0 # xt_pixel xt_is 54000 1.30 um spot 6
    #ronchi_sem_lib.ronchiC3Offset = -130.0 # xt_pixel xt_is 110000 1.5 um
    #ronchi_sem_lib.ronchiC3Offset = -30.0
    #sem.Pause('Please set C3 offset to where you can clearly see the global xLPP center')
    #ronchi_sem_lib.ronchiC3Offset = float(sem.ReportImageDistanceOffset()) - ronchi_sem_lib.ronchiStartC3Offset

def readCalibrations():
    global pixel_xt_matrix, is_xt_matrix, df_is_matrix
    # read calibrations from file. Only replace hardcoded default if
    # there is saved value
    r, ronchi_c3_offset = cal_util.readCalibration('pixel_xt_matrix', cal_dir, optics)
    if r:
        pixel_xt_matrix = r
    else:
        log(f'WARNING: no pixel_xt_matrix saved for {optics}, using default')
    r, ronchi_c3_offset = cal_util.readCalibration('is_xt_matrix', cal_dir, optics)
    if r:
        is_xt_matrix = r
    else:
        log(f'WARNING: no is_xt_matrix saved for {optics}, using default')
    r, dummy = cal_util.readCalibration('df_is_matrix', cal_dir)
    if r:
        df_is_matrix = r

def add_lpp_meta_to_next_mdoc():
    for k,v in (
            ('ImageDistanceOffset', sem.ReportImageDistanceOffset()),
        ):
        v_str = '%.12f' % (float(v))
        sem.AddToNextFrameStackMdoc(k, v_str)

def calc_xt_is(xt0, is_delta):
    # This form works for both array and list of list [x,y]
    xt1 = [0.0,0.0]
    xt1[0] = xt0[0]+is_delta[0]*is_xt_matrix[0][0]+is_delta[1]*is_xt_matrix[1][0]
    xt1[1] = xt0[1]+is_delta[0]*is_xt_matrix[0][1]+is_delta[1]*is_xt_matrix[1][1]
    return xt1

def calc_df_is(df0, is_delta):
    # This form works for both array and list of list
    df1 = 0.0
    df1 =df0+is_delta[0]*df_is_matrix[0][0]+is_delta[1]*df_is_matrix[1][1]
    return df1

def saveZeroImageShiftDefocusXLens():
    global lafisZeroImageShiftDefocus
    global lafisZeroImageShiftXLens
    lafisZeroImageShiftDefocus = sem.ReportDefocus()
    if ronchi_sem_lib.hasXLens:
        lafisZeroImageShiftXLens = [float(v) for v in sem.ReportXLensDeflector(2)[:2]]
    else:
        lafisZeroImageShiftXLens = None

def doLafis(is_x, is_y):
    global lafisIsDone, lafisXtCorrectionX, lafisXtCorrectionY
    log(f"WARNING: ***********doing LAFIS for image shift {is_x:.3f}, {is_y:.3f}")
    saveZeroImageShiftDefocusXLens()
    sem.AdjustBeamTiltforIS()
    df0 = lafisZeroImageShiftDefocus
    xt0 = lafisZeroImageShiftXLens
    is_delta = (is_x, is_y)
    df1 = calc_df_is(df0,is_delta)
    sem.SetDefocus(df1)
    if ronchi_sem_lib.hasXLens:
        xt1 = calc_xt_is(xt0,is_delta)
        sem.SetXLensDeflector(2, xt1[0], xt1[1])
        lafisXtCorrectionX = xt1[0] - xt0[0]
        lafisXtCorrectionY = xt1[1] - xt0[1]
    lafisIsDone = True

def restoreLafis():
    global lafisIsDone, lafisXtCorrectionX, lafisXtCorrectionY
    if not lafisIsDone:
        log("WARNING: LAFIS not done, can not restore")
        return
    sem.RestoreBeamTilt()
    sem.SetDefocus(lafisZeroImageShiftDefocus)
    if ronchi_sem_lib.hasXLens and lafisZeroImageShiftXLens is not None:
        xt_x, xt_y = lafisZeroImageShiftXLens
        sem.SetXLensDeflector(2, xt_x, xt_y)
        lafisXtCorrectionX = 0.0
        lafisXtCorrectionY = 0.0
    lafisIsDone = False
    log('WARNING: Lafis restored')

def _acquire_ronchi_image(trial_offset_baseline, ronchi_offset, pass_label=''):
    full_ronchi = ronchi_sem_lib.acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=pass_label)
    full_shape = full_ronchi.shape
    import mrcfile
    global count
    with mrcfile.new(f"ronchi{count:02d}.mrc", overwrite=True) as mrc:
        if np.issubdtype(full_ronchi.dtype, np.integer):
            mrc.set_data(full_ronchi.astype(np.float32))
        else:
            mrc.set_data(full_ronchi)
    count += 1
    #TODO corp off outside the beam if needed automatically
    return full_ronchi[:,int(0.25*full_shape[1]):]

def cross_correlate(img1, img2, shift=True):
    f1 = np.fft.fft2(img1)
    f2 = np.fft.fft2(img2)
    corr = np.fft.ifft2(f1 * np.conj(f2)).real
    if shift:
        corr = np.fft.fftshift(corr)  # swap quadrants: Q1<->Q3, Q2<->Q4
    return corr

def xlpp_center_finding(shifted_corr_arr, threshold_factor=0.5):
    """
    Find correlation peak by center of mass of a thresholded binary map.
    shifted_corr_arr should be a swapped/shifted correlation map where
    the center of the image is the origin if self correlated.
    Note: These are numpy convension with axes y,x
    """
    c_shape = shifted_corr_arr.shape
    c_center = c_shape[0]//2, c_shape[1]//2
    my_max = shifted_corr_arr.max()
    my_mean = shifted_corr_arr.mean()
    # center of mass in the thresholded correlation map gives a better
    # estimate of cc peak we should use.
    laser_threshold = my_mean+(my_max-my_mean)*threshold_factor
    laser=np.where(shifted_corr_arr > laser_threshold, 1, 0)
    laser_center = np.array(ndimage.center_of_mass(laser))
    # peak shift np array 
    corr_shift = np.array(c_center) - laser_center
    display_util.addImage(laser, [corr_shift,])
    return corr_shift

def _find_shift_numpy(img0, img1):
    """
    find shift of the highly off-plane xlpp image from img0 to img1
    return in numpy convention (row, col)
    """
    cor_image = cross_correlate(img0,img1, shift=True)
    # TODO need auto thresholding
    factor = 0.5 if is_simu else 0.8
    peak = xlpp_center_finding(cor_image, threshold_factor=factor)
    return peak

def _find_shift_sem(img0, img1):
    """
    find shift of the highly off-plane xlpp image from img0 to img1
    return in sem convention (col, row)
    """
    peak = _find_shift_numpy(img0, img1)
    peak_sem =  peak[::-1].copy()
    return peak_sem

def _measureLafisResidual(cal_image_shifts, trial_offset_baseline, ronchi_offset):
    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
    residual_shifts = np.zeros(cal_image_shifts.shape, dtype=np.float32)
    #display_util.addImage(img0_array)
    old_cal_xts = []
    for i, my_is in enumerate(cal_image_shifts):
        sem.SetImageShift(my_is[0],my_is[1])
        doLafis(my_is[0],my_is[1])
        img_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
        display_util.addImage(img_array)
        residual = _find_shift_sem(img0_array, img_array)
        residual_shifts[i] = np.array(residual)
        sem.SetImageShift(0,0)
        old_cal_xt_x = lafisXtCorrectionX
        old_cal_xt_y = lafisXtCorrectionY
        old_cal_xts.append([old_cal_xt_x,old_cal_xt_y])
        restoreLafis()
    print('lafis residual pixel shift x,y', residual_shifts)
    return np.array(old_cal_xts), residual_shifts

def update_is_xt_matrix(transform_arr):
    global is_xt_matrix
    is_xt_matrix = transform_arr.tolist()
    print('updated_is_xt_matrix', is_xt_matrix)

def _refineLafisMatrix(image_shift_scale, trial_offset_baseline, ronchi_offset):
    cal_image_shift_bases = np.array([[1,0],[0,1],[-1,0],[0,-1]])
    cal_image_shifts = image_shift_scale * cal_image_shift_bases
    print('image_shifts for cal', cal_image_shifts)
    old_cal_xts, pixel_residuals = _measureLafisResidual(cal_image_shifts, trial_offset_baseline, ronchi_offset)
    xt_residuals_arr = pixel_residuals @ np.array(pixel_xt_matrix)
    try:
        transform_matrix = cal_util.solveTransform(old_cal_xts-xt_residuals_arr, cal_image_shifts)
        update_is_xt_matrix(transform_matrix)
    except Exception as e:
        log(f'Error: Calibration not updated {e} Bad xt residuals {xt_residuals_arr}')
        transform_matrix = is_xt_matrix
    display_util.addTransformFitData(transform_matrix, old_cal_xts, xt_residuals_arr, cal_image_shifts)
    return pixel_residuals

def update_pixel_xt_matrix(transform_arr):
    if is_simu:
        # simulator correlation error is large. Don't update the value
        return
    global pixel_xt_matrix
    pixel_xt_matrix = transform_arr.tolist()

def _calibrate_pixel_xt_matrix(xt_scale, trial_offset_baseline, ronchi_c3_value):
    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    changes = np.array([[1.0,1.0],[1.0,-1.0],[-1,-1],[-1,1]])
    cal_xt_changes = xt_scale * changes
    pixel_shifts = changes * 0
    xt0 = np.array([float(v) for v in sem.ReportXLensDeflector(2)[:2]])
    for i, my_change in enumerate(cal_xt_changes):
        my_xt = my_change + xt0
        sem.SetXLensDeflector(2, my_xt[0], my_xt[1])
        log(f'xt set to {sem.ReportXLensDeflector(2)} in pixel_xt calibration')
        img_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
        my_shift = _find_shift_sem(img0_array, img_array)
        pixel_shifts[i] = np.array(my_shift)
        sem.SetXLensDeflector(2, xt0[0], xt0[1])
    try:
        observed_to_change_matrix = cal_util.solveTransform(cal_xt_changes, pixel_shifts)
        update_pixel_xt_matrix(observed_to_change_matrix)
    except Exception as e:
        log(f'Error: Calibration not updated {e} Bad pixel shift measured {pixel_shifts}')
        return None
    pixel_residuals = cal_xt_changes @ np.linalg.inv(observed_to_change_matrix) - pixel_shifts
    return pixel_residuals

def calibrateXtPixelMatrix():
    checkRonchigramSetup()
    saveZeroImageShiftDefocusXLens()
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    # Do pixel_xt_matrix calibration using an estimated xt_scale
    xt_scale = 0.000114648
    log(f'calibrating pixel_xt_matrix with xt change of {xt_scale} rad')
    print(f'calibrating pixel_xt_matrix with xt change of {xt_scale} rad')
    pixel_residuals = _calibrate_pixel_xt_matrix(xt_scale, trial_offset_baseline, ronchi_offset)
    if pixel_residuals is None:
        raise ValueError('pixel_xt_matrix calibration failed')
    if is_simu:
        # always accept
        return
    norm_0 = np.linalg.norm(pixel_xt_matrix[0])
    norm_1 = np.linalg.norm(pixel_xt_matrix[1])
    print(f'magnitudes: {norm_0:.3e}, {norm_1:.3e}')
    if abs(norm_0-norm_1) > 0.3 * (norm_0+norm_1)/2:
        raise ValueError('pixel_xt_matrix unreasonable. magnitudes between two axis too different')
    cos_theta = np.dot(pixel_xt_matrix[0], pixel_xt_matrix[1]) / (norm_0 * norm_1)
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    print(f'angle in degrees between the two axis {np.degrees(angle_rad):.2}')
    if abs(cos_theta) > np.sin(np.radians(10.0)):
        raise ValueError('pixel_xt_matrix unreasonable. angle between two axis too far from 90 degrees')
    print('pixel_residuals', pixel_residuals)

def calibrateLafis():
    checkRonchigramSetup()
    saveZeroImageShiftDefocusXLens()
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    cal_image_shift_scales = [1,2.5,5]    #in um
    max_trials = 3
    converging_deviation_threshold = 1.5e-5
    for i, scale in enumerate(cal_image_shift_scales):
        # do a refinement of the existing is_xt_matrix
        log(f'calibrating lafis_matrix with image shift of {scale} um')
        trial = 1
        while True:
            if trial > max_trials:
                raise ValueError('Lafis calibration did not converge.')
            log(f'calibrating lafis_matrix trial {trial:d} with image shift of {scale} um')
            pixel_residual = _refineLafisMatrix(scale, trial_offset_baseline, ronchi_offset)
            xt_residual = pixel_residual @ np.array(pixel_xt_matrix)
            mean_xt_deviation = np.linalg.norm(xt_residual).mean()
            if i == 0:
                # allow more tolerance on the first measurement
                this_threshold = converging_deviation_threshold * 2
            else:
                this_threshold = converging_deviation_threshold
            print(f'mean_deviation/scale={mean_xt_deviation/scale:.3e} vs converging_threshold={this_threshold:.3e}')
            if mean_xt_deviation /scale < this_threshold:

                break
            else:
                trial += 1
    resetOptics() 

def testXtPixel():
    xt0_x, xt0_y = lafisZeroImageShiftXLens
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    ronchi_c3_value = trial_offset_baseline + ronchi_offset
    cal_image_shift_scale = 0.5    #in um

    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    pixel_shifts = np.array([[0.0,0.0],[0.0,0.0]])
    quarter_x = img0_array.shape[0]//8
    pixel_shifts[1][1] = quarter_x
    xt0 = np.array((xt0_x, xt0_y))
    print('pixel_xt_matrix', pixel_xt_matrix)
    print('pixel_shift required',pixel_shifts)
    xt_delta = pixel_shifts @ np.array(pixel_xt_matrix)
    print('delta xt values will be applied',xt_delta)
    xt_total = xt0 + xt_delta
    sem.SetXLensDeflector(2,xt_total[1][0],xt_total[1][1])
    img1_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    resetOptics() 

def testLafis():
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = -20
    ronchi_c3_value = trial_offset_baseline + ronchi_offset
    all_is = [(0,0),(5,0),(-5,0),(0,5),(0,-5)]
    for my_is in all_is:
        sem.SetImageShift(my_is[0],my_is[1])
        doLafis(my_is[0],my_is[1])
        img1_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
        restoreLafis()
    resetOptics() 

if __name__=='__main__':
    #saveCalibrations()
    checkRonchigramSetup()
    print(ronchi_sem_lib.ronchiC3Offset)
    saveZeroImageShiftDefocusXLens()
    readCalibrations()
    ##### calibrate ronchiCorrMatrix
    failed_xt_pixel = False

    try:
        calibrateXtPixelMatrix()
        cal_util.saveCalibration('pixel_xt_matrix', cal_dir, session_name,pixel_xt_matrix, optics)
        if display_util.image_buffer:
            display_util.showImages()
            # reset image_buffer
            display_util.image_buffer = []
    except Exception as e:
        log(f'ERROR: Failed: {e}')
        log('Aborting....')
        failed_xt_pixel = True

    #testXtPixel()

    if not failed_xt_pixel:
        ##### lafis
        try:
            calibrateLafis()
            print(f'saving is_xt_matrix: {is_xt_matrix}')
            cal_util.saveCalibration('is_xt_matrix', cal_dir, session_name,is_xt_matrix, optics)
            if input('Perform test at 5 um IS shift ? (Y/N/y/n)').lower() == 'y':
                testLafis()
        except Exception as e:
            log(f'ERROR: Failed: {e}')
            print(f'Failed: {e}')
            log('Aborting....')
            print("""Tip:Try accurately align to on-plane""")
    print(f'final pixel_xt_matrix: {pixel_xt_matrix}')
    if not failed_xt_pixel:
        print(f'final is_xt_matrix: {is_xt_matrix}')
    if display_util.Fit_X:
        display_util.showTransformFit()
