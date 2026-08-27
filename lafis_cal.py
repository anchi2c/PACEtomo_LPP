#!Python
# ===================================================================
#ScriptName     Lafis Calibration
# Purpose:      Runs Calibration for LAFIS.
#               Make sure to run selectTargets script first to generate compatible Navigator settings and a target file.
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2021/04/16
# Revision:     v1.9.2c
# Last Change:  2026/05/27: selectable tilt schemes (dose_symmetric, bidirectional, continuous)
# ===================================================================

############ SETTINGS ############ 

########## Ronchigram / laser alignment ##########
# Trial LD area must match Record position; only exposure should differ.
# Overridable from target file via _bset (e.g. _bset doRonchigram true).
# Requires ronchi_sem_lib.hasXLens = True.
# Settings set in ronchi_sem_lib.py
########## END Ronchigram settings ##########

import sys
sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
import serialem as sem
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
xt_is_matrix = [[0.000324, -0.000347],[0.001100, 0.00028125]]  #26jul23
#xt_is_matrix = [[0.0001, -0.000347],[0.001100, 0.00028125]] #testing 
df_is_matrix = [[0.041381,0.012342], [0.041381,0.012342]]
xt_pixel_matrix = [[0.000100, -0.000],[0.000, 0.000100]]  #starting guess mrad/pixel

lafisZeroImageShiftDefocus = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisZeroImageShiftXLens = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisIsDone = False            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisXtCorrectionX = 0.0       # set from doLafis as the correction made on XLens
lafisXtCorrectionY = 0.0       # set from doLafis as the correction made on XLens
########## 

# ScriptName Script 11 Recall xt0 and other origin values from temp_xt0.json

def resetOptics():
    if platform.system() == 'Windows':
		# TODO: should use working directory
        filepath = 'X:\\k3f_serialem\\p25aug25a\\temp_xt0.json'
    else:
        filepath = './temp_xt0.json'
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
        color = 2
        style = 1 
    elif text.startswith("DEBUG:"):
        color = 1
        if breakpoints:
            breakpoint()
    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)

def saveCalibrations():
    #TODO save matrix as json
    cal_dir, session_name = cal_util.getCalibrationsDir()
    os.makedirs(cal_dir, exist_ok=True)
    mag,*_ = sem.ReportMag()
    cal_util.saveCalibration('xt_pixel_matrix_%d' % int(mag), cal_dir, session_name,xt_pixel_matrix) 
    cal_util.saveCalibration('xt_is_matrix', cal_dir, session_name,xt_is_matrix) 
    cal_util.saveCalibration('df_is_matrix', cal_dir, session_name,df_is_matrix) 

def readCalibrations():
    global xt_pixel_matrix, xt_is_matrix, df_is_matrix
    cal_dir, session_name = cal_util.getCalibrationsDir()
    mag,*_ = sem.ReportMag()
    # read calibrations from file. Only replace hardcoded default if
    # there is saved value
    r = cal_util.readCalibration('xt_pixel_matrix_%d' % int(mag), cal_dir)
    if r:
        xt_pixel_matrix = r
    r = cal_util.readCalibration('xt_is_matrix', cal_dir)
    if r:
        xt_is_matrix = r
    r = cal_util.readCalibration('df_is_matrix', cal_dir)
    if r:
        df_is_matrix = r

def add_lpp_meta_to_next_mdoc():
    for k,v in (
            ('ImageDistanceOffset', sem.ReportImageDistanceOffset()),
        ):
        v_str = '%.12f' % (float(v))
        sem.AddToNextFrameStackMdoc(k, v_str)

def checkRonchigramSetup():
    ronchi_sem_lib.checkRonchigramSetup()
    #ronchi_sem_lib.ronchiC3Offset = -173.0 # xt_pixel xt_is 88000 1.5 um
    ronchi_sem_lib.ronchiC3Offset = -130.0 # xt_pixel xt_is
    #ronchi_sem_lib.ronchiC3Offset = -30.0
    #sem.Pause('Please set C3 offset to where you can clearly see the global xLPP center')
    #ronchi_sem_lib.ronchiC3Offset = float(sem.ReportImageDistanceOffset()) - ronchi_sem_lib.ronchiStartC3Offset

def calc_xt_is(xt0, is_delta):
	# This form works for both array and list of list
    xt1 = [0.0,0.0]
    xt1[0] = xt0[0]+is_delta[0]*xt_is_matrix[0][0]+is_delta[1]*xt_is_matrix[1][0]
    xt1[1] = xt0[1]+is_delta[0]*xt_is_matrix[0][1]+is_delta[1]*xt_is_matrix[1][1]
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
        print('calculated xt',xt1)
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
    """
    shifted_corr_arr
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
    display_util.addImage(laser, corr_shift)
    return corr_shift

def find_shift(img0, img1):
    """
    find shift of the highly off-plane xlpp image from img0 to img1
    """
    cor_image = cross_correlate(img0,img1, shift=True)
    peak = xlpp_center_finding(cor_image, threshold_factor=0.8)
    print('shift on image', peak)
    return peak

def _measureLafisResidual(image_shift_scale, trial_offset_baseline, ronchi_offset):
    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
    residual_shifts = np.array([[0.0,0.0],[0.0,0.0]])
    cal_image_shifts = np.array([[image_shift_scale,0.0],[0.0,image_shift_scale]])
    cal_util.addImage(img0_array)
    for axis in (0,1):
        my_is = cal_image_shifts[axis]
        sem.SetImageShift(my_is[0],my_is[1])
        doLafis(my_is[0],my_is[1])
        img_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_offset)
        display_util.addImage(img_array)
        residual = find_shift(img0_array, img_array)
        residual_shifts[axis] = np.array(residual)
        sem.SetImageShift(0,0)
        restoreLafis()
    print('lafis residual', residual_shifts)
    return residual_shifts

def update_xt_is_matrix(xt_residual_arr, image_shift_scale):
    global xt_is_matrix
    xt_is_arr = np.array(xt_is_matrix) - xt_residual_arr / image_shift_scale
    xt_is_matrix = xt_is_arr.tolist()
    print('updated_xt_is_matrix', xt_is_matrix)

def _refineLafisMatrix(image_shift_scale, trial_offset_baseline, ronchi_offset):
    residuals = _measureLafisResidual(image_shift_scale, trial_offset_baseline, ronchi_offset)
    xt_residuals_arr = residuals @ np.linalg.inv(np.array(xt_pixel_matrix))
    #update_xt_is_matrix(xt_residuals_arr, image_shift_scale)
    return xt_residuals_arr

def update_xt_pixel_matrix(transform_arrx, shift_scale):
    global xt_pixel_matrix
    xt_is_arr = np.array(xt_pixel_matrix) + xt_residual_arr / shift_scale
    xt_is_matrix = xt_is_arr.tolist()

def _calibrate_xt_pixel_matrix(xt_scale, trial_offset_baseline, ronchi_c3_value):
    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    pixel_shifts = np.array([[0.0,0.0],[0.0,0.0]])
    cal_xt_shifts = np.array([[xt_scale,xt_scale],[-xt_scale,xt_scale]])
    xt0 = np.array([float(v) for v in sem.ReportXLensDeflector(2)[:2]])
    for axis in (0,1):
        my_xt = cal_xt_shifts[axis]
        sem.SetXLensDeflector(2, my_xt[0]+xt0[0], my_xt[1]+xt0[1])
        print(sem.ReportXLensDeflector(2))
        img_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
        my_shift = find_shift(img0_array, img_array)
        pixel_shifts[axis] = np.array(my_shift)
        sem.SetXLensDeflector(2, xt0[0], xt0[1])
    global xt_pixel_matrix
    try:
        xt_pixel_matrix = cal_util.solveTransform(cal_xt_shifts, pixel_shifts)
    except Exception as e:
        print('Error: Calibration not updated {e} Bad pixel shift measured {pixel_shifts}')
    return

def calibrateXtPixelMatrix():
    checkRonchigramSetup()
    saveZeroImageShiftDefocusXLens()
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    cal_image_shift_scale = 0.5    #in um

    # Step 1 do xt_pixel_matrix calibration using an estimated xt_scale
    # based on the xt_is_matrix and cal_image_shift_scale
    xt_scale = np.array(xt_is_matrix).mean() * cal_image_shift_scale
    print(f'calibrating xt_pixel_matrix with xt change of {xt_scale} rad')
    _calibrate_xt_pixel_matrix(xt_scale, trial_offset_baseline, ronchi_offset)

def calibrateLafis():
    checkRonchigramSetup()
    saveZeroImageShiftDefocusXLens()
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    cal_image_shift_scale = 2.5    #in um
    # Step 2 do a refinement of the existing xt_is_matrix
    print(f'calibrating lafis_matrix with image shift of {cal_image_shift_scale} um')
    _refineLafisMatrix(cal_image_shift_scale, trial_offset_baseline, ronchi_offset)
    resetOptics() 

def testXtPixel():
    xt0_x, xt0_y = lafisZeroImageShiftXLens
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
    ronchi_c3_value = trial_offset_baseline + ronchi_offset
    cal_image_shift_scale = 0.5    #in um

    img0_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    pixel_shifts = np.array([[0.0,0.0],[0.0,0.0]])
    quarter_x = img0_array.shape[1]//4
    pixel_shifts[1][1] = quarter_x
    xt0 = np.array((xt0_x, xt0_y))
    print('xt_pixel_matrix', xt_pixel_matrix)
    xt_delta = pixel_shifts @ np.array(xt_pixel_matrix)
    print('total xt values will be applied',xt0 + xt_delta)
    xt_total = xt0 + xt_delta
    sem.SetXLensDeflector(2,xt_total[1][0],xt_total[1][1])
    img1_array = _acquire_ronchi_image(trial_offset_baseline, ronchi_c3_value)
    resetOptics() 

def testLafis():
    trial_offset_baseline = ronchi_sem_lib.ronchiStartC3Offset
    ronchi_offset = ronchi_sem_lib.ronchiC3Offset
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
    saveZeroImageShiftDefocusXLens()
    #readCalibrations()
    # calibrate ronchiCorrMatrix
    calibrateXtPixelMatrix()
    testXtPixel()
    #calibrateLafis()  
    saveCalibrations()
    #testLafis()
    print(f'xt_pixel_matrix: {xt_pixel_matrix}')
    print(f'xt_is_matrix: {xt_is_matrix}')
    if display_util.image_buffer:
        display_util.showImages()
