#!Python
# ===================================================================
#ScriptName     PACEtomo_measureOffset
# Purpose:      Estimates tilt axis offset for PACEtomo. (Thanks to Wim Hagen for the suggestion!)
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2022/05/10
# Revision:     v1.8
# Last Change:  2026/05/27: CTF or beam-tilt measure; beam-tilt autofocus always (1.9)
# ===================================================================

############ SETTINGS ############ 

increment   = 5     # tilt step
maxTilt     = 15    # maximum +/- tilt angle
offset      = 5     # +/- offset for measured positions in microns from tilt axis (also accepts lists e.g. [2, 4, 6])

plot        = True  # plot measurements

# Defocus measurement in tilt loop (sem.G(-1) replacement)
defocusMethod = "ctf"       # ctf | beam_tilt

# X-tilt for CtfFind
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
ctfDefocusLo = -10.0        # CtfFind search range low [microns]
ctfDefocusHi = -0.2         # CtfFind search range high [microns]
ctf_resolution_max_A = 10.0 # retry CtfFind if resolution [A] is above this
ctf_max_attempts = 3        # max CtfFind attempts per measurement
ctf_retry_delay_s = 5       # delay before refocus shot on retry [s]

# Beam-tilt (match calibrate_beam_tilt_scaling.py)
tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
defocus_tilt_correction = beam_tilt_correction
beam_tilt_xtilt_x = 0.0
beam_tilt_xtilt_y = 0.0
spherical_aberration_mm = 2.7     # Cs for defocus = -disp/(2*beta) - Cs*beta^2 [mm]
autofocus_cycles = 2              # initial autofocus: correct focus to target_defocus (sem.G)
measure_cycles = 1                # cycles per tilt-loop measurement

target_defocus = -4.0             # defocus target for initial beam-tilt autofocus [microns]
target_defocus_tolerance_um = 0.05

########## END SETTINGS ########## 

import serialem as sem
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import PACEtomo_beamTiltDefocus as btdef

btdef.configure(sem_module=sem)

########### FUNCTIONS ###########

def dZ(alpha, y0):
    return y0 * np.tan(np.radians(-alpha))


def beam_tilt_measure_defocus():
    """Beam-tilt defocus via shared calibrated measurement."""
    return btdef.measure_defocus(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        defocus_tilt_correction=defocus_tilt_correction,
        xtilt_x=beam_tilt_xtilt_x,
        xtilt_y=beam_tilt_xtilt_y,
        cs_mm=spherical_aberration_mm,
    )


def ctf_measure_defocus():
    """CTF defocus: set X-tilt, Focus, CtfFind (with retries), restore X-tilt."""
    xtX, xtY = sem.ReportXLensDeflector(2)
    try:
        sem.SetXLensDeflector(2, ctfXtiltX, ctfXtiltY)
        cfind = []
        sem.NoMessageBoxOnError(1)
        try:
            for attempt in range(1, ctf_max_attempts + 1):
                if attempt > 1:
                    sem.Delay(ctf_retry_delay_s, "s")
                sem.F()
                cfind = sem.CtfFind("A", ctfDefocusLo, ctfDefocusHi)
                if len(cfind) == 0:
                    sem.Echo(f"ERROR: CtfFind failed on attempt {attempt}/{ctf_max_attempts}.")
                    if attempt < ctf_max_attempts:
                        continue
                    return np.nan
                resolution = float(cfind[-1])
                if resolution <= ctf_resolution_max_A:
                    break
                sem.Echo(
                    f"WARNING: CtfFind resolution {resolution:.2f} A > {ctf_resolution_max_A} A "
                    f"(attempt {attempt}/{ctf_max_attempts})"
                )
            else:
                sem.Echo(
                    f"WARNING: CtfFind resolution still > {ctf_resolution_max_A} A after "
                    f"{ctf_max_attempts} attempts; using last result."
                )
        finally:
            sem.NoMessageBoxOnError(0)
        defocus = float(cfind[0])
        sem.Echo(f"CtfFind: {defocus:.4f} microns ({cfind[-1]} A)")
        return defocus
    finally:
        sem.SetXLensDeflector(2, xtX, xtY)


def measure_defocus():
    """sem.G(-1): measure defocus only, no focus change."""
    if defocusMethod not in ("ctf", "beam_tilt"):
        sem.OKBox(f"ERROR: Unknown defocusMethod '{defocusMethod}'. Use 'ctf' or 'beam_tilt'.")
        sem.Exit()
    if defocusMethod == "beam_tilt":
        defocus = np.nan
        for _ in range(measure_cycles):
            defocus, speed_x, speed_y = beam_tilt_measure_defocus()
        sem.Echo(
            f"Defocus: {defocus:.4f} um, drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
        )
        return float(defocus)
    return float(ctf_measure_defocus())


def autofocus_apply():
    """sem.G: measure defocus and correct to target_defocus."""
    defocus = np.nan
    for cycle in range(1, autofocus_cycles + 1):
        defocus, speed_x, speed_y = beam_tilt_measure_defocus()
        error = target_defocus - defocus
        sem.Echo(
            f"Autofocus {cycle}/{autofocus_cycles}: "
            f"measured={defocus:.4f} um, target={target_defocus:.3f} um, error={error:.3f} um"
        )
        if abs(error) <= target_defocus_tolerance_um:
            return defocus
        sem.ChangeFocus(error)
    return defocus


def Tilt(tilt):
    sem.TiltTo(tilt)

    for i in range(len(offsets)):
        sem.ImageShiftByMicrons(0, offsets[i])
        sem.Delay(5, "s")
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
sem.Echo(f"Defocus measurement method (tilt loop): {defocusMethod}")
sem.Echo("Rough eucentricity...")
sem.GoToLowDoseArea("V")
sem.SetImageShift(0, 0)
sem.Eucentricity(1)
sem.GoToLowDoseArea("F")
sem.SetImageShift(0, 0)

sem.Echo(f"Beam-tilt autofocus to {target_defocus:.3f} microns (sem.G, always beam-tilt)...")
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