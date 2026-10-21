#!Python
import os
import numpy as np
import sys
import math
import platform
if platform.system() == 'Windows':
    is_simu = False
    sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
    import serialem as sem
else:
    is_simu = True
    print('testing on Mac/Linux with simulator')
    import sem_simulator as sem
import ronchi_cal
import ronchi_lib
import ronchi_sem_lib
import cal_util
import display_util

diff_ratio_limit = 0.005
angle_limit = 0.5

debug = False

ronchi_cal.ronchi_sem_lib.checkRonchigramSetup()
pixel_size_angs = ronchi_cal.sem.ReportCurrentPixelSize('T')
xt_tilt = 1 #(scaled at 1e-5 rad)
corr_scale = 1e-5
image_binning = 2
ronchi_binning = 32
peak_radius = 100
# The same value has always been used.
ronchi_fixed_pixel_size_um = 0.98e-4 * 2

x0_start, y0_start = sem.ReportCondenserStigmator()
x0 = x0_start
y0 = y0_start

delta = 0.001  # a known value to effectively change output as a step size.
n_step = 7     # this allows likely an interpolation of the good result.
steps = range((-n_step // 2)+1, (n_step // 2)+1)

def refineStigY(steps):
    """
    Stigmator Y mainly changes the angle between the 2 ronchigrams
    We want the angle to be 90 degrees.
    """
    stigys = []
    angles = []
    for i in steps:
        stigy = y0 + delta * i
        sem.SetCondenserStigmator(x0,stigy)
        pass_label = f'stigy change:{delta*i:.3e}'
        if is_simu:
            sem.test_rip = True # load -20 c3offset image for testing
        ks, phase = ronchi_cal.measure_ronchigram_ks_phases(ronchi_fixed_pixel_size_um, ronchi_binning,
                           peak_radius, corr_scale, pass_label)
        diff_ratio, angle_degree = ronchi_lib.analyze_ronchigram_stig(ks)
        if abs(diff_ratio) > 0.15 or abs(angle_degree-90) > 10:
            raise ValueError(f'bad {diff_ratio*100: .2f}%, {angle_degree}')
        stigys.append(stigy)
        angles.append(angle_degree)
    sem.SetCondenserStigmator(x0,y0)
    slopes, intercepts, residuals = cal_util.solveLines(np.array(stigys),np.array([angles,angles]).T)
    print(slopes, intercepts, residuals)
    right_angle_stigy = (90.0 - intercepts[0])/slopes[0]
    print(f'refined result for Y: {right_angle_stigy:.5f}')
    sem.SetCondenserStigmator(x0,right_angle_stigy)
    return right_angle_stigy


def refineStigX(steps):
    """
    Stigmator X mainly changes the relative magnitudes of the two ks.
    We want them to be equal length.
    """
    stigxs = []
    diffs = []
    for i in steps:
        stigx = x0 + delta * i
        sem.SetCondenserStigmator(stigx,y0)
        pass_label = f'stigx change:{delta*i:.3e}'
        if is_simu:
            sem.test_rip = True # load -20 c3offset image for testing
        ks, phase = ronchi_cal.measure_ronchigram_ks_phases(ronchi_fixed_pixel_size_um, ronchi_binning,
                           peak_radius, corr_scale, pass_label)
        diff_ratio, angle_degree = ronchi_lib.analyze_ronchigram_stig(ks)
        if abs(diff_ratio) > 0.15 or abs(angle_degree-90) > 10:
            raise ValueError(f'bad {diff_ratio*100: .2f}%, {angle_degree}')
        stigxs.append(stigx)
        diffs.append(diff_ratio)
    sem.SetCondenserStigmator(x0,y0)
    slopes, intercepts, residuals = cal_util.solveLines(np.array(stigxs),np.array([diffs,diffs]).T)
    print(slopes, intercepts, residuals)
    equal_magnitude_stigx = -intercepts[0]/slopes[0]
    print(f'refined result for X: {equal_magnitude_stigx:.5f}')
    sem.SetCondenserStigmator(equal_magnitude_stigx,y0)
    return equal_magnitude_stigx

def testStig(stigx,stigy):
    sem.SetCondenserStigmator(stigx,stigy)
    print('TESTING....')
    pass_label = f'stigx, stigy:{stigx:.5f}, {stigy:.5g}'
    if is_simu:
        sem.test_rip = True # load -20 c3offset image for testing
    ks, phase = ronchi_cal.measure_ronchigram_ks_phases(ronchi_fixed_pixel_size_um, ronchi_binning,
                           peak_radius, corr_scale, pass_label)
    diff_ratio, angle_degree = ronchi_lib.analyze_ronchigram_stig(ks)
    return abs(diff_ratio) < diff_ratio_limit and abs(angle_degree - 90.0) < angle_limit

trials = 1
max_trials = 3
while trials <= max_trials:
    if testStig(x0, y0):
        print('Alignment converged')
        break
    try:
        y0 = refineStigY(steps)
        print(f'current values: {x0:.5f}, {y0:.5f}')
        x0 = refineStigX(steps)
        print(f'current values: {x0:.5f}, {y0:.5f}')
        trials += 1
    except Exception as e:
        print(e)
        sem.SetCondenserStigmator(x0_start,y0_start)
        break
if trials > max_trials:
    print(f'Did not converge within {max_trials}')
    sem.SetCondenserStigmator(x0_start,y0_start)
