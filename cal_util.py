#!Python
# ===================================================================
#ScriptName     Calibration utilities
# Purpose:      Define and save calibrations in jsonl for persistent calibrations
# Author:       Anchi Cheng
# ===================================================================
import os
from datetime import datetime, timezone
import json
import numpy as np
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

timestampFormat = "%Y-%m-%d %H:%M:%S %Z"
illuminatedAreaTolerance = 1e-3   # in ReportIlluminatedArea unit
ronchiC3OffsetTolerance = 1.0   # in ReportImageDistanceOffset unit

def log(*args):
    print(args)
    pass

def getOpticsKey(ronchi_c3_offset=None):
    """
    Optics condition that calibrations such as pixel_xt_matrix depend on.
    ronchi_c3_offset is ronchi_sem_lib.ronchiC3Offset used to acquire
    the ronchigram.
    """
    mag,*_ = sem.ReportMag()
    spot_size = sem.ReportSpotSize()
    try:
        illuminated_area = float(sem.ReportIlluminatedArea())
    except Exception as e:
        # only available on some Thermo Scientific microscopes
        log('WARNING: illuminated area not available', e)
        illuminated_area = None
    if ronchi_c3_offset is not None:
        ronchi_c3_offset = float(ronchi_c3_offset)
    return {'mag': int(mag), 'spot_size': int(spot_size), 'illuminated_area': illuminated_area,
            'ronchi_c3_offset': ronchi_c3_offset}

def _valueMatch(v1, v2, tolerance):
    if v1 is None or v2 is None:
        return v1 is None and v2 is None
    return abs(v1 - v2) <= tolerance

def opticsMatch(optics1, optics2):
    if optics1['mag'] != optics2['mag'] or optics1['spot_size'] != optics2['spot_size']:
        return False
    if not _valueMatch(optics1['illuminated_area'], optics2['illuminated_area'], illuminatedAreaTolerance):
        return False
    # records saved before ronchi_c3_offset was added have no such key
    return _valueMatch(optics1.get('ronchi_c3_offset'), optics2.get('ronchi_c3_offset'), ronchiC3OffsetTolerance)

def saveCalibration(cal_type, cal_dir, session_name, data, optics=None):
    cal_path = os.path.join(cal_dir, cal_type+'.jsonl')
    log(cal_path)
    if isinstance(data, np.ndarray):
        data = data.tolist()
    cal_data = {}
    cal_data['timestamp'] = datetime.now().astimezone().strftime(timestampFormat)
    cal_data['session'] = session_name
    if optics is not None:
        cal_data['optics'] = optics
    cal_data['calibration'] = data
    # saved as JSONL: one JSON object per line
    with open(cal_path, "a") as f:
        f.write(json.dumps(cal_data)+"\n")
    return

def readCalibration(cal_type, cal_dir, optics=None):
    """
    Read the most recent calibration value from file. If optics is
    given, read the most recent one saved with matching optics.
    """
    cal_path = os.path.join(cal_dir, cal_type+'.jsonl')
    log('reading', cal_path)
    if not os.path.exists(cal_path):
        return None
    if optics is not None:
        return _readCalibrationMatchingOptics(cal_path, optics)
    # read from backward to get most recent entry
    with open(cal_path, "rb") as f:
        # go to the end position
        f.seek(0,2)
        pos = f.tell()
        line = b""
        while pos > 0:
            pos -= 1
            f.seek(pos)
            c = f.read(1)
            if c == b"\n" and line:
                break
            line = c + line
    my_data = json.loads(line.decode("utf-8"))
    if my_data:
        if 'timestamp' in my_data.keys():
            my_data['timestamp'] = datetime.strptime(my_data['timestamp'],timestampFormat)
        if 'calibration' in my_data.keys():
            return my_data['calibration']

def _readCalibrationMatchingOptics(cal_path, optics):
    """
    Search newest-first and stop at the first record with matching optics.
    Records saved without optics never match.
    """
    # text check to skip parsing lines of other mag
    mag_text = '"mag": %d,' % optics['mag']
    with open(cal_path, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if mag_text not in line:
            continue
        my_data = json.loads(line)
        if 'optics' in my_data and opticsMatch(my_data['optics'], optics):
            return my_data.get('calibration')
    return None

def getCalibrationsDir():
    working_dir = sem.ReportDirectory()
    root_dir, session_name = os.path.split(working_dir)
    if not session_name:
        root_dir, session_name = os.path.split(root_dir)

    cal_dir = os.path.join(root_dir,'calibrations')
    return cal_dir, session_name

def solveTransform(scope_changes, observed_shifts):
    """
    solve transformation array from observed_shifts to scope_changes
    """
    A = []
    B = []
    for (x,y), (xp, yp) in zip(scope_changes, observed_shifts):
        A.append([x,y, 0,0])
        A.append([0,0,x,y])
        B.append(xp)
        B.append(yp)
    data_src = np.array(A)
    data_results = np.array(B)
    try:
        params, residuals, rank, sv = np.linalg.lstsq(data_src, data_results, rcond=None)
        m11,m12,m21,m22 = params
        scope_to_observed = np.array([[m11,m12],[m21,m22]]).T
        return np.linalg.inv(scope_to_observed)
    except Exception as e:
        raise ValueError(f'Can not solve transform matrix. {e}')

def solveLines(x,y):
    """
    x is 1D array of n values
    y is 2D array of (n,m) values where m is differet dataset such
    as repeating measurement or independent axes.
    """
    print('solving x',x)
    print('against y',y)
    # matrix stays the same shape: (n_points, 2)
    A = np.vstack([x, np.ones(len(x))]).T
    # lstsq solves for all columns of y at once
    result, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)

    # result has shape [slope, intercept] pair per column of y
    slopes = result[0]      # array of m slopes
    intercepts = result[1]  # array of m intercepts
    return slopes, intercepts, residuals
