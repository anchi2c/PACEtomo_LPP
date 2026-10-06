#!Python
import os
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
import ronchi_sem_lib
import cal_util

def readRonchiCalibrations():
    # Get the latest calibration at the same optics except ronchiC3Offset
    optics = ronchi_sem_lib.getOpticsKey()
    ref_correct_ks, ronchi_c3_offset = cal_util.readCalibration('ronchi_ref_ks', cal_dir, optics)
    if ref_correct_ks is None:
        # No calibration to transfer.  Use module defaults
        return
    optics = ronchi_sem_lib.getOpticsKey(ronchi_c3_offset)  #limit to the same ronchi_c3_offset
    ref_phases, ronchi_c3_offset = cal_util.readCalibration('ronchi_ref_phase', cal_dir, optics)
    ronchi_sem_lib.ronchiC3Offset = ronchi_c3_offset
    ronchi_sem_lib.ronchiTargetPhaseA = ref_phases[0]           # vertical laser (rad)
    ronchi_sem_lib.ronchiTargetPhaseB = ref_phases[1]        # horizontal laser (rad)
    ronchi_sem_lib.ronchiCorrectKs    = ref_correct_ks

working_dir = sem.ReportDirectory()
cal_dir, session_name = cal_util.getCalibrationsDir(working_dir)
os.makedirs(cal_dir, exist_ok=True)

ronchi_sem_lib.checkRonchigramSetup()
readRonchiCalibrations()
print('phase A',ronchi_sem_lib.ronchiTargetPhaseA)
print('phase B',ronchi_sem_lib.ronchiTargetPhaseB)
print('correctKs', ronchi_sem_lib.ronchiCorrectKs)

ronchi_sem_lib.doRonchigramCorrection()
