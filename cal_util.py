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

def log(*args):
    print(args)
    pass

def saveCalibration(cal_type, cal_dir, session_name, data):
    cal_path = os.path.join(cal_dir, cal_type+'.jsonl')
    log(cal_path)
    if isinstance(data, np.ndarray):
        data = data.tolist()
    cal_data = {}
    cal_data['timestamp'] = datetime.now().astimezone().strftime(timestampFormat)
    cal_data['session'] = session_name
    cal_data['calibration'] = data
    # saved as JSONL: one JSON object per line
    with open(cal_path, "a") as f:
        f.write(json.dumps(cal_data)+"\n")
    return

def readCalibration(cal_type, cal_dir):
    """
    Read the most recent calibration value from file.
    """
    cal_path = os.path.join(cal_dir, cal_type+'.jsonl')
    log('reading', cal_path)
    if not os.path.exists(cal_path):
        return None
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
