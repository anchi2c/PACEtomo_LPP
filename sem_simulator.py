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

file_count = 0
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
cam_binning = 1
added_mdoc_dict = {}
test_failed = False

def pause():
    time.sleep(0.1)

def _log(name: str, *args, **kwargs) -> None:
    arg_str = ", ".join([repr(a) for a in args] + [f"{k}={v!r}" for k, v in kwargs.items()])
    #print(f"[sem stub] {name}({arg_str})")

def Echo(msg):
    _log(msg)
# ---------------------------------------------------------------------------
# Version / Session Control
# ---------------------------------------------------------------------------

def IsVersionAtLeast(*args, **kwargs) -> bool:
    _log("IsVersionAtLeast", *args, **kwargs)
    return True


def Exit(*args, **kwargs) -> None:
    _log("Exit", *args, **kwargs)


def IsVariableDefined(name: str) -> int:
    _log("IsVariableDefined", name)
    return 1


def SetPersistentVar(name: str, value: Any = "") -> None:
    _log("SetPersistentVar", name, value)


def GetVariable(name: str) -> str:
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
    return 1


def OKBox(message: str) -> None:
    _log("OKBox", message)


def Pause(message: str) -> None:
    _log("Pause", message)


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
    _log("ImageProperties", *args, **kwargs)
    return (0, 0, 1, 1.0)


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
    _log("CameraProperties", *args, **kwargs)
    return (0, 0)


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
    global file_count
    file_count += 1
    _log("R", *args, **kwargs)
    file_root = f"pt1_ts_{file_count:03d}"
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

def S(*args, **kwargs) -> None:
    _log("S", *args, **kwargs)
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)


def V(*args, **kwargs) -> None:
    _log("V", *args, **kwargs)


def L(*args, **kwargs) -> None:
    _log("L", *args, **kwargs)
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)


def F(*args, **kwargs) -> None:
    _log("F", *args, **kwargs)


def T(*args, **kwargs) -> None:
    _log("T", *args, **kwargs)
    global buffer_image_shift
    b_is_x, b_is_y, *_ = ReportImageShift() 
    buffer_image_shift = (b_is_x, b_is_y)


def AcquireToMatchBuffer(*args, **kwargs) -> None:
    _log("AcquireToMatchBuffer", *args, **kwargs)

import numpy as np
# ImageDistanceOffset at which the simulated fringes vanish (the value
# calibrate_ronchigram_start_c3 should recover).
ronchi_flat_c3_offset = -20.0
# ImageDistanceOffset at which the reference image below was acquired; it is
# reproduced unzoomed (zoom_factor == 1).
ref_c3_offset = -130.0
if not test_failed:
    sim_arr = mrc.read('./data/c3_offset_-130.mrc') #high underfocus for zoom and croping
    sq_size = int(min(sim_arr.shape)*0.8)
else:
    sim_arr = mrc.read('./data/lamella_failed.mrc') #ronchigram analysis would fail this one
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
    mag,*_ = ReportMag()
    cal_dir, session_name = cal_util.getCalibrationsDir()
    pixel_xt_matrix = cal_util.readCalibration('pixel_xt_matrix_%d' % int(mag), cal_dir)
    if pixel_xt_matrix is None:
        pixel_xt_matrix = np.array([[1.115e-6, 1.2858e-6],[5.903e-6, -3.3343e-7]])  #starting guess mrad/pixel
    else:
        pixel_xt_matrix = np.array(pixel_xt_matrix)
    pixel_is_matrix = np.array([[ 0.07404388, -0.02263578],
       [-0.01479404, -0.00302799]])
    #pixel_is_matrix = np.array([[ 0.06004388, -0.02263578],
    #   [-0.01479404, -0.00302799]])
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
        shift_x, shift_y = (shift_from_is + shift_from_xt).tolist()
        print('simu total pixel shift x,y', shift_x, shift_y)
        shift_arr = shift_pad_or_crop(sim_arr, int(shift_y), int(shift_x))
        # sem functions is X,Y
        cropped_arr = CropCenterToSize(shift_arr, cropped_shape[1], cropped_shape[0])
        c3_focus = ReportImageDistanceOffset()
        if not test_failed:
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


# ---------------------------------------------------------------------------
# Navigator
# ---------------------------------------------------------------------------

def ReportNavItem(*args, **kwargs) -> None:
    _log("ReportNavItem", *args, **kwargs)


def SetSelectedNavItem(*args, **kwargs) -> None:
    _log("SetSelectedNavItem", *args, **kwargs)


def MoveToNavItem(*args, **kwargs) -> None:
    _log("MoveToNavItem", *args, **kwargs)


def RealignToOtherItem(*args, **kwargs) -> None:
    _log("RealignToOtherItem", *args, **kwargs)


def LoadOtherMap(*args, **kwargs) -> None:
    _log("LoadOtherMap", *args, **kwargs)


def SaveNavigator(*args, **kwargs) -> None:
    _log("SaveNavigator", *args, **kwargs)


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


def SetFrameBaseName(*args, **kwargs) -> None:
    _log("SetFrameBaseName", *args, **kwargs)


def ReportFrameBaseName(*args, **kwargs) -> str:
    _log("ReportFrameBaseName", *args, **kwargs)
    return 0, 'test', 0


def SetFrameNameFormat(*args, **kwargs) -> None:
    _log("SetFrameNameFormat", *args, **kwargs)


def ReportFrameSavingPath(*args, **kwargs) -> str:
    _log("ReportFrameSavingPath", *args, **kwargs)
    return "NONE"


def ReportLastFrameFile(*args, **kwargs) -> tuple:
    _log("ReportLastFrameFile", *args, **kwargs)
    return ("", "", "")

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
    return 0


def SaveSettings(*args, **kwargs) -> None:
    _log("SaveSettings", *args, **kwargs)


def Delay(*args, **kwargs) -> None:
    _log("Delay", *args, **kwargs)


def LongOperation(*args, **kwargs) -> None:
    _log("LongOperation", *args, **kwargs)
