#!Python
# ===================================================================
#ScriptName     PACEtomo_measureOffset
# Purpose:      Estimates tilt axis offset for PACEtomo. (Thanks to Wim Hagen for the suggestion!)
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2022/05/10
# Revision:     v1.8
# Last Change:  2026/05/27: beam-tilt autofocus defocus measurement (1.9)
# ===================================================================

############ SETTINGS ############ 

increment   = 5     # tilt step
maxTilt     = 15    # maximum +/- tilt angle
offset      = 5     # +/- offset for measured positions in microns from tilt axis (also accepts lists e.g. [2, 4, 6])

plot        = True  # plot measurements

# Beam-tilt autofocus (calibrated correction factor)
tilt_angle_mrad = 5.0
beam_tilt_correction = 3 / 6.7
autofocus_cycles = 2   # cycles for initial autofocus (sem.G equivalent)
measure_cycles = 1     # cycles per defocus measurement (sem.G(-1) equivalent)

########## END SETTINGS ########## 

import serialem as sem
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

########### FUNCTIONS ###########

def dZ(alpha, y0):
    return y0 * np.tan(np.radians(-alpha))


def beam_tilt_measure_defocus():
    """
    Beam-tilt autofocus measurement (AutoFocus_New macro logic).
    Returns defocus [microns], drift speed x/y [nm/s].
    """
    beam_tilt = sem.ReportBeamTilt()
    tilt_x_orig = float(beam_tilt[0])
    tilt_y_orig = float(beam_tilt[1])
    tilt_x_plus = tilt_x_orig + beam_tilt_correction * tilt_angle_mrad
    tilt_x_minus = tilt_x_orig - beam_tilt_correction * tilt_angle_mrad

    pixel_size_binned = float(sem.ReportCurrentPixelSize("R"))
    binning = float(sem.ReportBinning("R"))
    pixel_size_unbinned = pixel_size_binned / binning

    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    sem.F()
    sem.ResetClock()
    sem.Copy("A", "L")

    sem.SetBeamTilt(tilt_x_minus, tilt_y_orig)
    sem.F()
    sem.AlignTo("L", 1)
    align_shift_1 = sem.ReportAlignShift()
    disp_x1_px = float(align_shift_1[0])
    disp_y1_px = float(align_shift_1[1])

    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    sem.F()
    elapsed = float(sem.ReportClock())

    sem.SetBeamTilt(tilt_x_orig, tilt_y_orig)
    sem.AlignTo("L", 1)
    align_shift_2 = sem.ReportAlignShift()
    disp_x2_px = float(align_shift_2[0])
    disp_y2_px = float(align_shift_2[1])

    drift_x = disp_x2_px * pixel_size_unbinned
    drift_y = disp_y2_px * pixel_size_unbinned
    speed_x = drift_x / elapsed if elapsed > 0 else 0.0
    speed_y = drift_y / elapsed if elapsed > 0 else 0.0

    displacement_from_tilt_x = (disp_x1_px - disp_x2_px / 2.0) * pixel_size_unbinned
    displacement_from_tilt_y = (disp_y1_px - disp_y2_px / 2.0) * pixel_size_unbinned
    displacement = np.sqrt(
        displacement_from_tilt_x * displacement_from_tilt_x
        + displacement_from_tilt_y * displacement_from_tilt_y
    )

    if displacement_from_tilt_x == 0:
        sign = 1.0
    else:
        sign = displacement_from_tilt_x / abs(displacement_from_tilt_x)

    defocus_measured = -1.0 * sign * displacement / tilt_angle_mrad
    return defocus_measured, speed_x, speed_y


def measure_defocus():
    """Measure defocus only (sem.G(-1) equivalent)."""
    sem.L()
    defocus = np.nan
    for _ in range(measure_cycles):
        defocus, speed_x, speed_y = beam_tilt_measure_defocus()
    sem.Echo(
        f"Beam-tilt defocus: {defocus:.4f} microns, "
        f"drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
    )
    return float(defocus)


def autofocus_apply():
    """Measure defocus and apply focus correction (sem.G equivalent)."""
    defocus = np.nan
    for _ in range(autofocus_cycles):
        defocus, speed_x, speed_y = beam_tilt_measure_defocus()
        sem.ChangeFocus(-defocus)
    sem.Echo(
        f"Beam-tilt autofocus applied, last defocus: {defocus:.4f} microns, "
        f"drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
    )


def Tilt(tilt):
    sem.TiltTo(tilt)

    for i in range(len(offsets)):
        sem.ImageShiftByMicrons(0, offsets[i])
        focus[i].append(measure_defocus())
        sem.SetImageShift(0, 0)

    if tilt == 0:
        for j in range(len(offsets)): 
            focus0.append(focus[j][-1])

    angles.append(float(tilt))
    
###########################

sem.ResetClock()
sem.SuppressReports()
sem.SetUserSetting("ShiftToTiltAxis", 1, 1)

oldOffset = sem.ReportTiltAxisOffset()[0]

sem.Echo("Currently set tilt axis offset: " + str(oldOffset))

sem.Echo("##### Starting tilt axis offset estimation #####")
sem.Echo("Rough eucentricity...")
sem.Eucentricity(1)

sem.Echo("Beam-tilt autofocus...")
autofocus_apply()

sem.Echo("Start tilt series...")
starttilt = -maxTilt
sem.TiltTo(starttilt)
sem.TiltBy(-increment)

offsets = [0]
if isinstance(offset, (list, tuple)):
    for val in offset:
        offsets.extend([-val, val])
else:
    offsets.extend([-offset, offset])

angles = []
focus = [[] for i in range(len(offsets))]
focus0 = []

steps = 2 * maxTilt / increment + 1

tilt = starttilt
for i in range(int(steps)):
    sem.Echo("Tilt to " + str(tilt) + " deg")
    Tilt(tilt)
    tilt += increment

relFocus = focus
for i in range(len(angles)):
    for j in range(len(offsets)):
        relFocus[j][i] -= focus0[j]

y0 = np.zeros(len(offsets))
y0_neg = np.zeros(len(offsets))
y0_pos = np.zeros(len(offsets))
for j in range(len(offsets)):
    y0[j], cov = optimize.curve_fit(dZ, angles, relFocus[j], p0=0)
    y0_neg[j], cov = optimize.curve_fit(dZ, [angle for angle in angles if angle <= 0], relFocus[j][:len([angle for angle in angles if angle <= 0])], p0=0)
    y0_pos[j], cov = optimize.curve_fit(dZ, [angle for angle in angles if angle >= 0], relFocus[j][len([angle for angle in angles if angle < 0]):], p0=0)

sem.Echo("Remaining tilt axis offsets:")
for i in range(0, len(offsets)):
    sem.Echo("[" + str(offsets[i]) + "]: " + str(round(y0[i] + offsets[i], 2)) + " (neg: " + str(round(y0_neg[i] + offsets[i], 2)) + ", pos: " + str(round(y0_pos[i] + offsets[i], 2)) + ")")
avgOffset = sum(y0) / len(offsets)
avgOffset_neg = sum(y0_neg) / len(offsets)
avgOffset_pos = sum(y0_pos) / len(offsets)
sem.Echo("Average remaining tilt axis offset: " + str(round(avgOffset, 2)) + " (neg: " + str(round(avgOffset_neg, 2)) + ", pos: " + str(round(avgOffset_pos, 2)) + ")")
totalOffset = round(avgOffset + oldOffset, 2)
sem.Echo("##############################################")
sem.Echo("Estimated total tilt axis offset: " + str(totalOffset))
sem.Echo("##############################################")

sem.TiltTo(0)
sem.ResetImageShift()

sem.SuppressReports(0)
sem.ReportClock()

if plot:
    offsets, relFocus = zip(*sorted(zip(offsets, relFocus)))    # ensure right order for plot points
    fig = plt.figure(figsize=(8, 6), tight_layout=True)
    plt.title('Z Shifts [microns]')
    for i in range(len(angles)):
        values = []
        for j in range(len(offsets)):
            values.append(relFocus[j][i])
        plt.plot(offsets, values, label=str(angles[i]) + " deg")

    plt.legend()
    plt.show()

userInput = sem.YesNoBox("The estimated total tilt axis offset is " + str(totalOffset) + ". Do you want to set the new tilt axis offset?")
if userInput == 1:
    sem.SetTiltAxisOffset(totalOffset)
    sem.Echo("The new tilt axis offset has been set!")
sem.Exit()