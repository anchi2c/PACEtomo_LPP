#!/usr/bin/env python
"""
Stub / mock implementation of the SerialEM Python scripting interface
(the `serialem` module normally provided by SerialEM at runtime).

This module is intended for local development, linting, and testing of
scripts written against SerialEM's `sem` API when SerialEM itself is not
running. Every method is a stub: it prints the call it received and
returns a plausible placeholder value based on how it is used in
practice. Replace with the real SerialEM `serialem` module when running
against an actual microscope.

Writen with help by Claude
"""

from typing import Any
import time
from pyami import mrc
import numpy as np
import scipy.ndimage as ndimage
import os
import sys

spot_size = 1
illuminated_area = 0.013023 # 1 um is 0.01 in sem
c3_offset = 0.0
x_tilt = (0.0,0.0)
image_shift = (0.0,0.0)
beam_tilt = (0.0,0.0)
defocus = 0.0
buffer_image_shift = (0.0,0.0)
buffer_x_tilt = (0.0,0.0)
stage_tilt = 0.0
stage_xyz = (0.0,0.0,0.0)
mag = 88000
pixel_size = 1.0  # um
cam_binning = 2
frame_saving_path = "./"
frame_base_name = ""
sequential_frame_num = 1
frame_name_format_angle = False
frame_name_format_sequential_num = True

added_mdoc_dict = {}
# script variables set by SetVariable and by sem calls such as ReportNavItem
sim_variables = {}


## when test_rip is True, a -20 c3offset file is loaded
test_rip = False

def pause():
    time.sleep(0.1)

def _log(name: str, *args, **kwargs) -> None:
    arg_str = ", ".join([repr(a) for a in args] + [f"{k}={v!r}" for k, v in kwargs.items()])
    print(f"[sem stub] {name}({arg_str})")

def Echo(msg):
    _log(msg)
# ---------------------------------------------------------------------------
# Version / Session Control
# ---------------------------------------------------------------------------

def IsVersionAtLeast(*args, **kwargs) -> bool:
    #_log("IsVersionAtLeast", *args, **kwargs)
    return True


def Exit(*args, **kwargs) -> None:
    _log("Exit", *args, **kwargs)


def IsVariableDefined(name: str) -> int:
    _log("IsVariableDefined", name)
    return 1


def SetPersistentVar(name: str, value: Any = "") -> None:
    _log("SetPersistentVar", name, value)


def SetVariable(name: str, value: Any = "") -> None:
    _log("SetVariable", name, value)
    sim_variables[name] = value


def listToSEMarray(values) -> str:
    """python list to SerialEM array variable string"""
    return "\n".join(str(v) for v in values)


def GetVariable(name: str) -> str:
    if name in sim_variables:
        return sim_variables[name]
    if name == 'navIndex':
        return 1
    if name == 'navNote':
        return '/Users/anchi.cheng/packages/sem-scripts/pt1_tgts.txt'
    _log("GetVariable", name)
    return ""


def ResetClock(*args, **kwargs) -> None:
    _log("ResetClock", *args, **kwargs)


def ReportClock(*args, **kwargs) -> float:
    _log("ReportClock", *args, **kwargs)
    return 0.0


def SuppressReports(*args, **kwargs) -> None:
    _log("SuppressReports", *args, **kwargs)


def ProgramTimeStamps(*args, **kwargs) -> None:
    _log("ProgramTimeStamps", *args, **kwargs)


# ---------------------------------------------------------------------------
# Dialogs / User Interaction
# ---------------------------------------------------------------------------

def YesNoBox(message: str) -> int:
    _log("YesNoBox", message)
    answer = input(message).lower() == 'y'
    return int(answer)


def OKBox(message: str) -> None:
    _log("OKBox", message)


def Pause(message: str) -> None:
    _log("Pause", message)


def EnterDefaultedNumber(*args, **kwargs) -> float:
    """default, number of decimal places, prompt; returns the default"""
    _log("EnterDefaultedNumber", *args, **kwargs)
    return args[0]


def EnterString(*args, **kwargs) -> None:
    """variable name, prompt; sets the variable to 'sim'"""
    _log("EnterString", *args, **kwargs)
    sim_variables[args[0]] = 'sim'


# ---------------------------------------------------------------------------
# Column / Gun / FEG
# ---------------------------------------------------------------------------

def AreDewarsFilling(*args, **kwargs) -> int:
    _log("AreDewarsFilling", *args, **kwargs)
    return 0


def IsFEGFlashingAdvised(*args, **kwargs) -> int:
    _log("IsFEGFlashingAdvised", *args, **kwargs)
    return 0


def NextFEGFlashHighTemp(*args, **kwargs) -> None:
    _log("NextFEGFlashHighTemp", *args, **kwargs)


def ReportColumnOrGunValve(*args, **kwargs) -> int:
    _log("ReportColumnOrGunValve", *args, **kwargs)
    return 1


def SetColumnOrGunValve(*args, **kwargs) -> None:
    _log("SetColumnOrGunValve", *args, **kwargs)


def ReportSpotSize(*args, **kwargs) -> int:
    _log("ReportSpotSize", *args, **kwargs)
    return spot_size

def SetSpotSize(*args, **kwargs) -> None:
    _log("SetSpotSize", *args, **kwargs)
    global spot_size
    spot_size = args[0]


def ReportIlluminatedArea(*args, **kwargs) -> float:
    _log("ReportIlluminatedArea", *args, **kwargs)
    return illuminated_area

def SetIlluminatedArea(*args, **kwargs) -> None:
    _log("SetIlluminatedArea", *args, **kwargs)
    global illuminated_area
    illuminated_area = args[0]


# ---------------------------------------------------------------------------
# Stage / Tilt
# ---------------------------------------------------------------------------

def ReportTiltAngle(*args, **kwargs) -> float:
    _log("ReportTiltAngle", *args, **kwargs)
    return stage_tilt


def TiltTo(*args, **kwargs) -> None:
    _log("TiltTo", *args, **kwargs)
    global stage_tilt
    stage_tilt = args[0]
    pause()

def TiltBy(*args, **kwargs) -> None:
    _log("TiltBy", *args, **kwargs)
    global stage_tilt
    stage_tilt += args[0]
    pause()

def ReportStageXYZ(*args, **kwargs) -> tuple:
    _log("ReportStageXYZ", *args, **kwargs)
    return stage_xyz


def MoveStageTo(*args, **kwargs) -> None:
    _log("MoveStageTo", *args, **kwargs)
    global stage_xyz
    if len(args) == 3:
        stage_xyz = (args[0],args[1],args[2])
    elif len(args) == 2:
        stage_xyz = (args[0],args[1],stage_xyz[2])

def ReportTiltAxisOffset(*args, **kwargs) -> tuple:
    _log("ReportTiltAxisOffset", *args, **kwargs)
    return (0.0, 0.0)


def Eucentricity(*args, **kwargs) -> None:
    _log("Eucentricity", *args, **kwargs)


def UpdateItemZ(*args, **kwargs) -> None:
    _log("UpdateItemZ", *args, **kwargs)


# ---------------------------------------------------------------------------
# Image Shift / Specimen Shift
# ---------------------------------------------------------------------------

def SetImageShift(*args, **kwargs) -> None:
    _log("SetImageShift", *args, **kwargs)
    global image_shift
    image_shift = (args[0], args[1])


def GetImageShift(*args, **kwargs) -> tuple:
    _log("GetImageShift", *args, **kwargs)
    return image_shift


def ReportImageShift(*args, **kwargs) -> tuple:
    _log("ReportImageShift", *args, **kwargs)
    print('reported image shift',image_shift)
    return image_shift


def ImageShiftByMicrons(*args, **kwargs) -> None:
    _log("ImageShiftByMicrons", *args, **kwargs)
    global image_shift
    image_shift = (image_shift[0]+args[0], image_shift[1]+args[1])

def ResetImageShift(*args, **kwargs) -> None:
    """move stage by the image shift and set image shift to zero"""
    _log("ResetImageShift", *args, **kwargs)
    global image_shift, stage_xyz
    # IS to specimen and specimen to stage are identity in the simulator
    stage_xyz = (stage_xyz[0] + image_shift[0], stage_xyz[1] + image_shift[1], stage_xyz[2])
    image_shift = (0.0, 0.0)


def ImageShiftByUnits(*args, **kwargs) -> None:
    _log("ImageShiftByUnits", *args, **kwargs)


def ImageShiftByPixels(*args, **kwargs) -> None:
    _log("ImageShiftByPixels", *args, **kwargs)
    global image_shift
    pixel_size = 0.0032
    image_shift = (image_shift[0]+args[0]*pixel_size, image_shift[1]+args[1]*pixel_size)


def ReportSpecimenShift(*args, **kwargs) -> tuple:
    _log("ReportSpecimenShift", *args, **kwargs)
    return (0.0, 0.0)


def ReportISforBufferShift(*args, **kwargs) -> tuple:
    _log("ReportISforBufferShift", *args, **kwargs)
    return buffer_image_shift


def ReportBeamTilt(*args, **kwargs) -> tuple:
    _log("ReportBeamTilt", *args, **kwargs)
    return beam_tilt

def SetBeamTilt(*args, **kwargs) -> None:
    _log("SetBeamTilt", *args, **kwargs)
    global beam_tilt
    beam_tile = args[0],args[1]

def ReportObjectiveStigmator(*args, **kwargs) -> tuple:
    _log("ReportObjectiveStigmator", *args, **kwargs)
    return (0.0, 0.0)

def SetObjectiveStigmator(*args, **kwargs) -> None:
    _log("SetObjectiveStigmator", *args, **kwargs)


def LimitNextAutoAlign(*args, **kwargs) -> None:
    _log("LimitNextAutoAlign", *args, **kwargs)


def SetImageDistanceOffset(*args, **kwargs) -> None:
    _log("SetImageDistanceOffset", *args, **kwargs)
    global c3_offset
    c3_offset = args[0] 

def ReportImageDistanceOffset(*args, **kwargs) -> float:
    _log("ReportImageDistanceOffset", *args, **kwargs)
    return c3_offset


# ---------------------------------------------------------------------------
# Focus / Defocus / Beam Tilt / CTF
# ---------------------------------------------------------------------------


def SetDefocus(*args, **kwargs) -> None:
    _log("SetDefocus", *args, **kwargs)
    global defocus
    defocus = args[0]

def ReportDefocus(*args, **kwargs) -> float:
    _log("ReportDefocus", *args, **kwargs)
    return defocus


def ChangeFocus(*args, **kwargs) -> None:
    _log("ChangeFocus", *args, **kwargs)
    global defocus
    defocus = defocus + args[0]


def SetTargetDefocus(*args, **kwargs) -> None:
    _log("SetTargetDefocus", *args, **kwargs)


def RefineZLP(*args, **kwargs) -> None:
    _log("RefineZLP", *args, **kwargs)


def AdjustBeamTiltforIS(*args, **kwargs) -> None:
    _log("AdjustBeamTiltforIS", *args, **kwargs)
    global beam_tilt
    beam_tilt = (beam_tilt[0] + image_shift[0]*0.1, beam_tilt[1] + image_shift[1]*0.1)

def RestoreBeamTilt(*args, **kwargs) -> None:
    _log("RestoreBeamTilt", *args, **kwargs)
    global beam_tilt
    beam_tilt = (0.0,0.0)

def ReportXLensDeflector(*args, **kwargs) -> tuple:
    _log("ReportXLensDeflector", *args, **kwargs)
    return x_tilt


def SetXLensDeflector(*args, **kwargs) -> None:
    _log("SetXLensDeflector", *args, **kwargs)
    global x_tilt
    x_tilt = (args[1],args[2])
    pause()

def ReportLDDefocusOffset(*args, **kwargs) -> float:
    _log("ReportLDDefocusOffset", *args, **kwargs)
    return 0.0


def CtfFind(*args, **kwargs) -> tuple:
    _log("CtfFind", *args, **kwargs)
    return (0.0,)


def Ctfplotter(*args, **kwargs) -> tuple:
    _log("Ctfplotter", *args, **kwargs)
    return (0.0,)


def ReportComaVsISmatrix(*args, **kwargs) -> tuple:
    _log("ReportComaVsISmatrix", *args, **kwargs)
    return (0.1, 0.0, 0.0, 0.1)


# ---------------------------------------------------------------------------
# Autofocus / Autotuning
# ---------------------------------------------------------------------------

def G(*args, **kwargs)-> None:
    # autofocus
    _log("G", *args, **kwargs)

def ReportAutoFocus(*args, **kwargs) -> list:
    _log("ReportDefocus", *args, **kwargs)
    return (ReportDefocus(), 0)


# ---------------------------------------------------------------------------
# Low Dose / Camera Areas
# ---------------------------------------------------------------------------

def GoToLowDoseArea(*args, **kwargs) -> None:
    _log("GoToLowDoseArea", *args, **kwargs)


def UpdateLowDoseParams(*args, **kwargs) -> None:
    _log("UpdateLowDoseParams", *args, **kwargs)


def RestoreLowDoseParams(*args, **kwargs) -> None:
    _log("RestoreLowDoseParams", *args, **kwargs)


def SetCameraArea(*args, **kwargs) -> None:
    _log("SetCameraArea", *args, **kwargs)


def RestoreCameraSet(*args, **kwargs) -> None:
    _log("RestoreCameraSet", *args, **kwargs)


def ReportLDAreaShift(*args, **kwargs) -> tuple:
    _log("ReportLDAreaShift", *args, **kwargs)
    return (0.0, 0.0)


def ReportAxisPosition(*args, **kwargs) -> tuple:
    _log("ReportAxisPosition", *args, **kwargs)
    return (0,)


# ---------------------------------------------------------------------------
# Exposure / Camera Properties
# ---------------------------------------------------------------------------

def ReportExposure(*args, **kwargs) -> tuple:
    _log("ReportExposure", *args, **kwargs)
    return (1.0, 0.0)


def SetExposure(*args, **kwargs) -> None:
    _log("SetExposure", *args, **kwargs)


def SetBinning(*args, **kwargs) -> None:
    _log("SetBinning", *args, **kwargs)
    global cam_binning
    cam_binning = args[0]


def ImageProperties(*args, **kwargs) -> tuple:
    """width, height, binning, exposure, pixel size (nm), mode"""
    _log("ImageProperties", *args, **kwargs)
    return (sq_size, sq_size, cam_binning, 1.0, pixel_size * 1000, 1)


def ImageConditions(*args, **kwargs) -> tuple:
    _log("ImageConditions", *args, **kwargs)
    return (0.0,)


def ReportMeanCounts(*args, **kwargs) -> float:
    _log("ReportMeanCounts", *args, **kwargs)
    return 0.0


def ReportCurrentPixelSize(*args, **kwargs) -> float:
    _log("ReportCurrentPixelSize", *args, **kwargs)
    # this is a preset camera binned pixel in micron
    preset = args[0]
    return pixel_size


def ReportCameraProperty(*args, **kwargs) -> float:
    _log("ReportCameraProperty", *args, **kwargs)
    return 0.0


def CameraProperties(*args, **kwargs) -> tuple:
    """unbinned width, height, rotation/flip, physical pixel size (um), pixel size at current mag (nm)"""
    _log("CameraProperties", *args, **kwargs)
    return (sq_size * cam_binning, sq_size * cam_binning, 0, 5.0, pixel_size * 1000 / cam_binning)


def ReportBinning(*args, **kwargs) -> int:
    _log("CameraProperties", *args, **kwargs)
    return cam_binning


def SetProperty(*args, **kwargs) -> None:
    _log("SetProperty", *args, **kwargs)


def ReportProperty(*args, **kwargs) -> float:
    _log("ReportProperty", *args, **kwargs)
    name = args[0]
    if name == 'MaximumTiltAngle':
        return 70
    return 0.0


def SetUserSetting(*args, **kwargs) -> None:
    _log("SetUserSetting", *args, **kwargs)


# ---------------------------------------------------------------------------
# Mag
# ---------------------------------------------------------------------------

def ReportMag(*args, **kwargs) -> tuple:
    _log("ReportMag", *args, **kwargs)
    return (mag,)


def SetMag(*args, **kwargs) -> None:
    _log("SetMag", *args, **kwargs)
    global mag
    mag = args[0]

# ---------------------------------------------------------------------------
# Acquisition (Trial / Record / View / Focus / Preview)
# ---------------------------------------------------------------------------

def R(*args, **kwargs) -> None:
    global sequential_frame_num
    _log("R", *args, **kwargs)
    file_root = f"{ReportFrameBaseName()[1]}"
    if frame_name_format_sequential_num:
        file_root += f"_{sequential_frame_num:03d}"
    if frame_name_format_angle:
        file_root += f"_{ReportTiltAngle():03.2f}"
    f=open(file_root+".mrc.mdoc",'w')
    pre_text = '''
PixelSpacing = 3.2
Voltage = 300
Version = SerialEM Version 4.2.24 64-bit,  built Jun 14 2026  10:00:20
ImageFile = pt717_ts_001.mrc
ImageSize = 2880 2048
DataMode = 1
'''
    keys = list(added_mdoc_dict.keys())
    keys.sort()
    added_list = list(map((lambda x: f'{x} = {added_mdoc_dict[x]}'), keys))
    my_text = '\n'.join(added_list) + '\n'
    post_text = ''' 
[T = SerialEM: CZIi Krios 2 XL-on                            13-Aug-26  12:45:07]

[T = Tilt axis angle = 156.5, binning = 2  spot = 6  camera = 1]

[ZValue = 0]
MinMaxMean = 30 2174 826.578
TiltAngle = -0.0124786
StagePosition = 467.928 473.526
StageZ = 80.2568
'''
    text = pre_text + my_text +  post_text
    f.write(text)
    f.close()
    mrc.write(np.zeros((4,4)),file_root+'.mrc')
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)
    #Save fake tif file
    frame_dir = ReportFrameSavingPath()
    base_name = f"{ReportFrameBaseName()[1]}"
    if frame_name_format_sequential_num:
        base_name += f"_{sequential_frame_num:03d}"
        sequential_frame_num += 1
    if frame_name_format_angle:
        base_name += f"_{ReportTiltAngle():03.2f}"
    file_name = base_name+".tif"
    full_path = os.path.join(frame_dir,file_name)
    print(full_path)
    if os.path.exists(full_path):
        print(f'Existing path: {full_path}')
        sys.exit()
    f=open(full_path,'w')
    f.write('fake tif\n')
    f.close()

def S(*args, **kwargs) -> None:
    _log("S", *args, **kwargs)
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)


def V(*args, **kwargs) -> None:
    _log("V", *args, **kwargs)


def Search(*args, **kwargs) -> None:
    _log("Search", *args, **kwargs)


def L(*args, **kwargs) -> None:
    _log("L", *args, **kwargs)
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)


def F(*args, **kwargs) -> None:
    _log("F", *args, **kwargs)


def T(*args, **kwargs) -> None:
    _log("T", *args, **kwargs)
    global sequential_frame_num, buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)
    frame_dir = ReportFrameSavingPath()
    base_name = f"{ReportFrameBaseName()[1]}"
    if frame_name_format_sequential_num:
        base_name += f"_{sequential_frame_num:03d}"
        sequential_frame_num += 1
    if frame_name_format_angle:
        base_name += f"_{ReportTiltAngle():03.2f}"
    file_name = base_name+".tif"
    full_path = os.path.join(frame_dir,file_name)
    print(full_path)
    if os.path.exists(full_path):
        print(f'Existing path: {full_path}')
        sys.exit()
    f=open(full_path,'w')
    f.write('fake tif\n')
    f.close()


def AcquireToMatchBuffer(*args, **kwargs) -> None:
    _log("AcquireToMatchBuffer", *args, **kwargs)

import numpy as np
# ImageDistanceOffset at which the simulated fringes vanish (the value
# calibrate_ronchigram_start_c3 should recover).
ronchi_flat_c3_offset = 0.0
# xt correction per um of defocus (the value defocus_xt_cal.py should recover)
sim_df_xt_vector = np.array([2.0e-4, -1.0e-4])
if not test_rip:
    # ImageDistanceOffset at which the reference image below was acquired; it is
    # reproduced unzoomed (zoom_factor == 1).
    ref_c3_offset = -130.0
    sim_arr = mrc.read('./data/c3_offset_-130.mrc') #high underfocus for zoom and croping
    sq_size = int(min(sim_arr.shape)*0.8)
else:
    ref_c3_offset = -20.0
    sim_arr = mrc.read('./data/c3_offset_-20.mrc') #ronchigram analysis should give [[9,848, -0.779],[0.779,9.381]]
    sq_size = int(min(sim_arr.shape))
cropped_shape = (sq_size, sq_size)

def shift_pad_or_crop(img, pad_y, pad_x):
    # from Claude with modification so the output shape is the same shape as the input
    # handle y-axis
    img_shape = img.shape
    if pad_y >= 0:
        img = np.pad(img, ((pad_y, 0), (0, 0)))
        img = img[0:img_shape[0],:]
    else:
        img = img[-pad_y:, :]
        img = np.pad(img, ((0, -pad_y), (0, 0)))
    # handle x-axis
    if pad_x >= 0:
        img = np.pad(img, ((0, 0), (pad_x, 0)))
        img = img[:,0:img_shape[1]]
    else:
        img = img[:, -pad_x:]
        img = np.pad(img, ((0, 0), (0,-pad_x)))
    return img

def zoom_center_to_size(data, zoom_factor=1):
    ny, nx = data.shape

    # size of the crop window, shrinks as zoom_factor grows
    crop_h = int(ny / zoom_factor)
    crop_w = int(nx / zoom_factor)
    cy, cx = ny // 2, nx // 2
    y0 = max(cy - crop_h // 2, 0)
    y1 = min(y0 + crop_h, ny)
    x0 = max(cx - crop_w // 2, 0)
    x1 = min(x0 + crop_w, nx)

    cropped = data[y0:y1, x0:x1]
    # interpolate back up to original resolution for a smooth zoomed view
    order = 1
    scale_y = ny / cropped.shape[0]
    scale_x = nx / cropped.shape[1]
    zoomed = ndimage.zoom(cropped, (scale_y, scale_x), order=order)

    return zoomed

def bufferImage(*args, **kwargs) -> Any:
    _log("bufferImage", *args, **kwargs)
    """
    Retrieve image array. Tom simulate image shift and x tilt, we report
    their values and scale the value to create a shift.
    """
    import cal_util
    import ronchi_sem_lib
    working_dir = ReportDirectory()
    cal_dir, session_name = cal_util.getCalibrationsDir(working_dir)
    pixel_xt_matrix, ronchi_c3_offset = cal_util.readCalibration('pixel_xt_matrix', cal_dir, ronchi_sem_lib.getOpticsKey(ronchi_c3_offset=ronchi_sem_lib.ronchiC3Offset))
    if pixel_xt_matrix is None:
        pixel_xt_matrix, ronchi_c3_offset = cal_util.readCalibration('pixel_xt_matrix', cal_dir, ronchi_sem_lib.getOpticsKey(ronchi_c3_offset=None))
        if pixel_xt_matrix is None:
            pixel_xt_matrix = np.array([[1.338e-6, -2.006312e-7],[2.36091e-7, 1.5751567e-6]])  #starting guess mrad/pixel
    else:
        pixel_xt_matrix = np.array(pixel_xt_matrix)
    print('pixel_xt_matrix', pixel_xt_matrix)

    # calculate pixel_is_matrix from saved is_xt_matrix
    is_xt_matrix, ronchi_c3_offset = cal_util.readCalibration('is_xt_matrix', cal_dir, ronchi_sem_lib.getOpticsKey(ronchi_c3_offset=None))
    pixel_is_matrix = -np.linalg.inv(is_xt_matrix) @ pixel_xt_matrix
    print('perfect pixel_is_matrix', pixel_is_matrix)
    # for lafis high off-plane
    is_scale = 1
    xt_scale = 1
    shift_x = 0
    shift_y = 0
    try:
        is_x,is_y,*_ = ReportImageShift()
        xt_x,xt_y,*_ = ReportXLensDeflector(2)
        shift_from_xt = np.array((xt_x,xt_y)) @ np.linalg.inv(pixel_xt_matrix) * xt_scale
        shift_from_is = np.array((is_x,is_y)) @ np.linalg.inv(pixel_is_matrix) * is_scale
        print('shift_from_xt', xt_x, xt_y, shift_from_xt)
        print('shift_from_is', is_x, is_y, shift_from_is)
        # defocus shifts xlpp the same way as xt, opposite to its correction
        shift_from_df = -defocus * sim_df_xt_vector @ np.linalg.inv(pixel_xt_matrix) * xt_scale
        shift_x, shift_y = (shift_from_is + shift_from_xt + shift_from_df).tolist()
        print('simu total pixel shift x,y', shift_x, shift_y)
        shift_arr = shift_pad_or_crop(sim_arr, int(shift_y), int(shift_x))
        # sem functions is X,Y
        cropped_arr = CropCenterToSize(shift_arr, cropped_shape[1], cropped_shape[0])
        c3_focus = ReportImageDistanceOffset()
        if not test_rip:
            # distance from the fringe-free offset; zoom is the ratio of the
            # reference image distance to the current one.
            c3_distance = c3_focus - ronchi_flat_c3_offset
            # avoid division 0
            if c3_distance == 0:
                c3_distance = 0.2
            zoom_factor = abs((ref_c3_offset - ronchi_flat_c3_offset)/c3_distance)
        else:
            zoom_factor = 1
        print('****zoom_factor',zoom_factor)
        zoomed_arr = zoom_center_to_size(cropped_arr, zoom_factor)
        return zoomed_arr
    except Exception as e:
        _log(f"Exception {e}")
        return [[0]]


def CropCenterToSize(*args, **kwargs) -> None:
    #_log("CropCenterToSize", *args, **kwargs)
    _log("CropCenterToSize", args[1], args[2])
    x = args[1] #x dimension
    y = args[2] #y dimension
    my_image = args[0]
    my_shape = my_image.shape
    center = my_shape[0]//2,my_shape[1]//2
    if x > my_shape[1] or y > my_shape[0]:
        raise ValueError(f'SEM CropCenterToSize can not crop img of shape (x,y)= ({my_shape[1]:d},{my_shape[0]:d}) to ({x:d},{y:d})')
    crop_start_y = max(center[0]-y//2,0)
    crop_start_x = max(center[1]-x//2,0)
    crop_end_y = crop_start_y + y
    crop_end_x = crop_start_x + x
    return my_image[crop_start_y:crop_end_y,crop_start_x:crop_end_x]

# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def AlignTo(*args, **kwargs) -> None:
    _log("AlignTo", *args, **kwargs)


def AlignBetweenMags(*args, **kwargs) -> None:
    _log("AlignBetweenMags", *args, **kwargs)


def ReportAlignShift(*args, **kwargs) -> tuple:
    _log("ReportAlignShift", *args, **kwargs)
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def Copy(*args, **kwargs) -> None:
    _log("Copy", *args, **kwargs)


def AddBufToStackWindow(*args, **kwargs) -> None:
    _log("AddBufToStackWindow", *args, **kwargs)


def ReportCurrentBuffer(*args, **kwargs) -> str:
    _log("ReportCurrentBuffer", *args, **kwargs)
    return "A"


def ReduceImage(*args, **kwargs) -> None:
    """buffer, binning; result goes to buffer A. Image is not reduced in the simulator"""
    _log("ReduceImage", *args, **kwargs)


def SnapshotToFile(*args, **kwargs) -> None:
    _log("SnapshotToFile", *args, **kwargs)


# ---------------------------------------------------------------------------
# Navigator
# ---------------------------------------------------------------------------

# Simulated navigator table. Indices used by the sem API are 1-based.
# Item 1 matches the navNote returned by GetVariable before any nav call.
nav_items = [
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'type': 0, 'group': 0, 'acquire': 0, 'color': 0, 'draw': 1,
     'label': '1', 'note': '/Users/anchi.cheng/packages/sem-scripts/pt1_tgts.txt'},
]
nav_current_index = 1
nav_next_group_id = 1000
nav_file = os.path.abspath('./sim_navigator.nav')
# nav item type: 0 point, 1 polygon, 2 map
NAV_POINT, NAV_POLYGON, NAV_MAP = 0, 1, 2

def _nav_index(index):
    """Convert sem index (1-based, 0 = current, negative = from end) to list index."""
    index = int(index)
    if index == 0:
        index = nav_current_index
    elif index < 0:
        index = len(nav_items) + 1 + index
    if not 1 <= index <= len(nav_items):
        raise IndexError(f'[sem stub] Navigator index {index} out of range 1-{len(nav_items)}')
    return index - 1

def _nav_add(x, y, z, item_type, group=0, note=''):
    nav_items.append({'x': float(x), 'y': float(y), 'z': float(z), 'type': item_type,
                      'group': int(group), 'acquire': 0, 'color': 0, 'draw': 1,
                      'label': str(len(nav_items) + 1), 'note': note})
    return len(nav_items)

def _nav_report(i):
    """Set nav variables as SerialEM does and return (index, x, y, z, type)."""
    item = nav_items[i]
    sim_variables['navLabel'] = item['label']
    sim_variables['navNote'] = item['note']
    sim_variables['navIndex'] = i + 1
    sim_variables['navAcqIndex'] = item['acquire']
    return (i + 1, item['x'], item['y'], item['z'], item['type'])


def ReportNavItem(*args, **kwargs) -> tuple:
    _log("ReportNavItem", *args, **kwargs)
    return _nav_report(_nav_index(0))


def ReportOtherItem(*args, **kwargs) -> tuple:
    _log("ReportOtherItem", *args, **kwargs)
    return _nav_report(_nav_index(args[0]))


def ReportNumTableItems(*args, **kwargs) -> int:
    _log("ReportNumTableItems", *args, **kwargs)
    return len(nav_items)


def ReportGroupStatus(*args, **kwargs) -> tuple:
    """index of current item, its group ID, number of items in the group"""
    _log("ReportGroupStatus", *args, **kwargs)
    i = _nav_index(0)
    group = nav_items[i]['group']
    n_in_group = len([item for item in nav_items if group and item['group'] == group])
    return (i + 1, group, n_in_group)


def ReportIfNavOpen(*args, **kwargs) -> int:
    """0 no navigator, 1 open, 2 open with a file"""
    _log("ReportIfNavOpen", *args, **kwargs)
    return 2


def ReportNavFile(*args, **kwargs) -> str:
    _log("ReportNavFile", *args, **kwargs)
    return nav_file


def ReportItemAcquire(*args, **kwargs) -> int:
    _log("ReportItemAcquire", *args, **kwargs)
    return nav_items[_nav_index(args[0] if args else 0)]['acquire']


def SetItemAcquire(*args, **kwargs) -> None:
    _log("SetItemAcquire", *args, **kwargs)
    value = int(args[1]) if len(args) > 1 else 1
    nav_items[_nav_index(args[0] if args else 0)]['acquire'] = value


def ChangeItemNote(*args, **kwargs) -> None:
    _log("ChangeItemNote", *args, **kwargs)
    nav_items[_nav_index(args[0])]['note'] = str(args[1]) if len(args) > 1 else ''


def ChangeItemLabel(*args, **kwargs) -> None:
    _log("ChangeItemLabel", *args, **kwargs)
    nav_items[_nav_index(args[0])]['label'] = str(args[1])


def ChangeItemColor(*args, **kwargs) -> None:
    _log("ChangeItemColor", *args, **kwargs)
    nav_items[_nav_index(args[0])]['color'] = int(args[1])


def ChangeItemDraw(*args, **kwargs) -> None:
    """toggle drawing if no value is given"""
    _log("ChangeItemDraw", *args, **kwargs)
    item = nav_items[_nav_index(args[0])]
    item['draw'] = int(args[1]) if len(args) > 1 else 1 - item['draw']


def DeleteNavigatorItem(*args, **kwargs) -> None:
    _log("DeleteNavigatorItem", *args, **kwargs)
    global nav_current_index
    i = _nav_index(args[0])
    del nav_items[i]
    if nav_current_index > len(nav_items):
        nav_current_index = max(1, len(nav_items))


def NavIndexWithNote(*args, **kwargs) -> int:
    """index of the first item with the note, 0 if none"""
    _log("NavIndexWithNote", *args, **kwargs)
    for i, item in enumerate(nav_items):
        if item['note'] == args[0]:
            return i + 1
    return 0


def NavIndexWithLabel(*args, **kwargs) -> int:
    """index of the first item with the label, 0 if none"""
    _log("NavIndexWithLabel", *args, **kwargs)
    for i, item in enumerate(nav_items):
        if item['label'] == str(args[0]):
            return i + 1
    return 0


def GetUniqueNavID(*args, **kwargs) -> int:
    _log("GetUniqueNavID", *args, **kwargs)
    global nav_next_group_id
    nav_next_group_id += 1
    return nav_next_group_id


def AddStagePosAsNavPoint(*args, **kwargs) -> int:
    """x, y, z [group ID]; returns index of the new point"""
    _log("AddStagePosAsNavPoint", *args, **kwargs)
    group = args[3] if len(args) > 3 else 0
    return _nav_add(args[0], args[1], args[2], NAV_POINT, group)


def AddImagePosAsNavPoint(*args, **kwargs) -> int:
    """buffer, image x, y [z] [group ID]; returns index of the new point"""
    _log("AddImagePosAsNavPoint", *args, **kwargs)
    width, height, _, _, pixel_nm, *_ = ImageProperties(args[0])
    stage_x, stage_y, stage_z = stage_xyz
    # camera to stage is identity in the simulator (see CameraToSpecimenMatrix)
    x = stage_x + (float(args[1]) - width / 2) * pixel_nm / 1000
    y = stage_y + (float(args[2]) - height / 2) * pixel_nm / 1000
    group = args[4] if len(args) > 4 else 0
    return _nav_add(x, y, stage_z, NAV_POINT, group)


def AddStagePointsAsPolygon(*args, **kwargs) -> int:
    """names of variables holding x and y arrays, z; returns index of the new polygon"""
    _log("AddStagePointsAsPolygon", *args, **kwargs)
    xs = [float(v) for v in str(sim_variables[args[0]]).split()]
    ys = [float(v) for v in str(sim_variables[args[1]]).split()]
    index = _nav_add(np.mean(xs), np.mean(ys), args[2], NAV_POLYGON)
    nav_items[index - 1]['pts'] = (xs, ys)
    return index


def GetNavGroupStageCoords(*args, **kwargs) -> None:
    """group ID, names of variables to set with x, y, z arrays of the group items"""
    _log("GetNavGroupStageCoords", *args, **kwargs)
    group_items = [item for item in nav_items if item['group'] == int(args[0])]
    for key, var_name in zip(('x', 'y', 'z'), args[1:4]):
        sim_variables[var_name] = listToSEMarray([item[key] for item in group_items])


def NewMap(*args, **kwargs) -> int:
    """[0] [note]; map at current stage position, returns its index"""
    _log("NewMap", *args, **kwargs)
    note = str(args[1]) if len(args) > 1 else ''
    return _nav_add(*stage_xyz, NAV_MAP, note=note)


def AdjustStagePosForNav(*args, **kwargs) -> tuple:
    """stage x, y, IS x, IS y -> stage position of the image-shifted area"""
    _log("AdjustStagePosForNav", *args, **kwargs)
    # IS to specimen and specimen to stage are identity in the simulator
    return (float(args[0]) + float(args[2]), float(args[1]) + float(args[3]))


def SetSelectedNavItem(*args, **kwargs) -> None:
    _log("SetSelectedNavItem", *args, **kwargs)
    global nav_current_index
    nav_current_index = _nav_index(args[0]) + 1


def MoveToNavItem(*args, **kwargs) -> None:
    _log("MoveToNavItem", *args, **kwargs)


def RealignToOtherItem(*args, **kwargs) -> None:
    _log("RealignToOtherItem", *args, **kwargs)


def LoadOtherMap(*args, **kwargs) -> None:
    _log("LoadOtherMap", *args, **kwargs)


def SaveNavigator(*args, **kwargs) -> None:
    """write nav_items as a minimal navigator autodoc file"""
    _log("SaveNavigator", *args, **kwargs)
    global nav_file
    if args:
        nav_file = os.path.abspath(args[0])
    lines = ['AdocVersion = 2.00', '']
    for item in nav_items:
        lines += [f"[Item = {item['label']}]",
                  f"Color = {item['color']}",
                  f"StageXYZ = {item['x']} {item['y']} {item['z']}",
                  f"Type = {item['type']}",
                  f"GroupID = {item['group']}",
                  f"Acquire = {item['acquire']}",
                  f"Draw = {item['draw']}"]
        if item['note']:
            lines.append(f"Note = {item['note']}")
        if 'pts' in item:
            xs, ys = item['pts']
            lines += [f"NumPts = {len(xs)}",
                      "PtsX = " + " ".join(str(v) for v in xs),
                      "PtsY = " + " ".join(str(v) for v in ys)]
        lines.append('')
    with open(nav_file, 'w') as f:
        f.write("\n".join(lines))


# ---------------------------------------------------------------------------
# Files / Frames
# ---------------------------------------------------------------------------

def OpenOldFile(*args, **kwargs) -> None:
    _log("OpenOldFile", *args, **kwargs)


def OpenNewFile(*args, **kwargs) -> None:
    _log("OpenNewFile", *args, **kwargs)


def CloseFile(*args, **kwargs) -> None:
    _log("CloseFile", *args, **kwargs)


def ReportFileNumber(*args, **kwargs) -> int:
    _log("ReportFileNumber", *args, **kwargs)
    return 0


def ReadFile(*args, **kwargs) -> None:
    _log("ReadFile", *args, **kwargs)


def ReadOtherFile(*args, **kwargs) -> None:
    _log("ReadOtherFile", *args, **kwargs)


def ReportFileZsize(*args, **kwargs) -> int:
    _log("ReportFileZsize", *args, **kwargs)
    return 1


def SetNewFileType(*args, **kwargs) -> None:
    _log("SetNewFileType", *args, **kwargs)


def AllowFileOverwrite(*args, **kwargs) -> None:
    _log("AllowFileOverwrite", *args, **kwargs)


def DoesFileExist(*args, **kwargs) -> int:
    _log("DoesFileExist", *args, **kwargs)
    return int(os.path.exists(args[0]))


def ReportFrameBaseName(*args, **kwargs) -> str:
    _log("ReportFrameBaseName", *args, **kwargs)
    return 1, frame_base_name, 0

def SetFrameBaseName(*args, **kwargs) -> None:
    _log("SetFrameBaseName", *args, **kwargs)
    global frame_base_name, sequential_frame_num
    frame_base_name = args[3]
    sequential_frame_num = 1


def SetFrameNameFormat(*args, **kwargs) -> None:
    _log("SetFrameNameFormat", *args, **kwargs)
    global frame_name_format_num, frame_name_format_angle
    if args[2] == 0x40:
        # turn off Sequential number
        frame_name_format_sequential_num = False
    if args[2] == 0x400:
        # turn on tilt angle
        frame_name_format_angle = True


def SetFrameSavingPath(*args, **kwargs) -> None:
    _log("SetFrameSavingPath", *args, **kwargs)
    global frame_saving_path
    frame_saving_path = args[0]

def ReportFrameSavingPath(*args, **kwargs) -> str:
    _log("ReportFrameSavingPath", *args, **kwargs)
    return frame_saving_path


def SetSequentialFrameNum(*args, **kwargs) -> None:
    _log("SetFrameSavingPath", *args, **kwargs)
    global sequential_frame_num
    sequential_frame_num = args[0]

def ReportSequentialFrameNum(*args, **kwargs) -> str:
    _log("ReportFrameSavingPath", *args, **kwargs)
    return sequential_frame_num


def ReportLastFrameFile(*args, **kwargs) -> tuple:
    _log("ReportLastFrameFile", *args, **kwargs)
    dirname = ReportFrameSavingPath()
    basename = ReportFrameBaseName()[1]
    fullname = basename+f'_{ReportSequentialFrameNum():03d}.tif'
    return (os.path.join(dirname,fullname), dirname, basename)

def AddToNextFrameStackMdoc(*args, **kwargs) -> None:
    _log("AddToNextFrameStackMdoc", *args, **kwargs)
    added_mdoc_dict[args[0]] = args[1]

# ---------------------------------------------------------------------------
# Coordinate Matrices
# ---------------------------------------------------------------------------

def StageToSpecimenMatrix(*args, **kwargs) -> tuple:
    _log("StageToSpecimenMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


def ISToSpecimenMatrix(*args, **kwargs) -> tuple:
    _log("ISToSpecimenMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


def SpecimenToISMatrix(*args, **kwargs) -> tuple:
    _log("SpecimenToISMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


def CameraToSpecimenMatrix(*args, **kwargs) -> tuple:
    _log("CameraToSpecimenMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


def SpecimenToStageMatrix(*args, **kwargs) -> tuple:
    _log("SpecimenToStageMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


def SpecimenToCameraMatrix(*args, **kwargs) -> tuple:
    _log("SpecimenToCameraMatrix", *args, **kwargs)
    return (1.0, 0.0, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Dose Tracking
# ---------------------------------------------------------------------------

def AreaForCumulRecordDose(*args, **kwargs) -> None:
    _log("AreaForCumulRecordDose", *args, **kwargs)


def AccumulateRecordDose(*args, **kwargs) -> None:
    _log("AccumulateRecordDose", *args, **kwargs)


# ---------------------------------------------------------------------------
# Autodoc / Log
# ---------------------------------------------------------------------------

def AddToAutodoc(*args, **kwargs) -> None:
    _log("AddToAutodoc", *args, **kwargs)


def WriteAutodoc(*args, **kwargs) -> None:
    _log("WriteAutodoc", *args, **kwargs)


def SaveLog(*args, **kwargs) -> None:
    _log("SaveLog", *args, **kwargs)


def SaveLogOpenNew(*args, **kwargs) -> None:
    _log("SaveLogOpenNew", *args, **kwargs)


def SetNextLogOutputStyle(*args, **kwargs) -> None:
    _log("SetNextLogOutputStyle", *args, **kwargs)


def EchoBreakLines(*args, **kwargs) -> None:
    _log("EchoBreakLines", *args, **kwargs)
    print('Script Log:', args[0])

def SetStatusLine(*args, **kwargs) -> None:
    _log("SetStatusLine", *args, **kwargs)


def ClearStatusLine(*args, **kwargs) -> None:
    _log("ClearStatusLine", *args, **kwargs)


# ---------------------------------------------------------------------------
# Directory
# ---------------------------------------------------------------------------

def ReportDirectory(*args, **kwargs) -> str:
    _log("ReportDirectory", *args, **kwargs)
    return "./"


def UserSetDirectory(*args, **kwargs) -> None:
    _log("UserSetDirectory", *args, **kwargs)


# ---------------------------------------------------------------------------
# Misc / Error Handling
# ---------------------------------------------------------------------------

def NoMessageBoxOnError(*args, **kwargs) -> None:
    _log("NoMessageBoxOnError", *args, **kwargs)


def KeyBreak(*args, **kwargs) -> int:
    _log("KeyBreak", *args, **kwargs)
    return 1


def SaveSettings(*args, **kwargs) -> None:
    _log("SaveSettings", *args, **kwargs)


def Delay(*args, **kwargs) -> None:
    _log("Delay", *args, **kwargs)


def LongOperation(*args, **kwargs) -> None:
    _log("LongOperation", *args, **kwargs)
