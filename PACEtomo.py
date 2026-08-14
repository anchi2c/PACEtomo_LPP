#!Python
# ===================================================================
#ScriptName     PACEtomo
# Purpose:      Runs parallel dose-symmetric tilt series on many targets with geometrical predictions using the first target as tracking tilt series.
#               Make sure to run selectTargets script first to generate compatible Navigator settings and a target file.
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2021/04/16
# Revision:     v1.9.2c
# Last Change:  2026/05/27: selectable tilt schemes (dose_symmetric, bidirectional, continuous)
# ===================================================================

############ SETTINGS ############ 

startTilt       = 0         # starting tilt angle [degrees] (should be divisible by step)
minTilt         = -45       # minimum absolute tilt angle [degrees]
maxTilt         = 45        # maximum absolute tilt angle [degrees]
step            = 3         # tilt step [degrees]
groupSize       = 2         # group size for dose_symmetric scheme (contiguous images per side before switching)
tiltScheme      = "dose_symmetric"  # dose_symmetric | bidirectional | continuous
minDefocus      = -0.5        # minimum defocus [microns] of target range (low defocus)
maxDefocus      = -1.5        # maximum defocus [microns] of target range (high defocus)
stepDefocus     = 0.5       # step [microns] between target defoci (between TS)

focusSlope      = 0.0       # [DEPRECATED] empirical linear focus correction [microns per degree] (obtained by linear regression of CTF fitted defoci over tilt series; microscope stage dependent)
delayIS         = 2.0      # delay [s] between applying image shift and Record
delayTilt       = 2.0       # delay [s] after stage tilt
zeroExpTime     = 0         # set to exposure time [s] used for start tilt image, if 0: use same exposure time for all tilt images
zeroDefocus     = 0         # set to defocus [microns] used for start tilt image, if 0: use same defocus for all tilt images

nav_item_list   = [1]        # e.g. [5, 10, 15]; empty = current nav item only (SetSelectedNavItem, SerialEM 4.2+)
nav_pretilt_list = []       # parallel to nav_item_list; fall back to global pretilt when empty
nav_rotation_list = []      # parallel to nav_item_list; fall back to global rotation when empty
nav_start_defocus_list = [] # parallel to nav_item_list; objective defocus [um] before autofocus; else default_nav_start_defocus
default_nav_start_defocus = -4  # SetDefocus at start of each nav item when not listed

# Track settings
trackExpTime    = 0         # set to exposure time [s] used for tracking tilt series, if 0: use same exposure time for all tilt series
trackDefocus    = 0         # set to defocus [microns] used for tracking tilt series, if 0: use same defocus range for all tilt series
trackMag        = 0         # set to nominal magnification for tracking tilt series (make sure detector is still covered under the same beam conditions), if 0: use same mag for all tilt series
trackTwice      = False     # track in 2 steps, useful when large tracking shifts cause inaccuracies in alignment and hence, in residual errors for all targets, but causes double exposure of tracking area

# Geometry settings
pretilt         = 0         # pretilt [degrees] of sample in deg e.g. after FIB milling (if milling direction is not perpendicular to the tilt axis, estimate and add rotation)
rotation        = 0         # rotation [degrees] of lamella vs tilt axis in deg (should be 0 deg if lamella is oriented perpendicular to tilt axis)
measureGeo      = False     # estimates pretilt and rotation values of sample by measuring defocus on geo points saved in target setup or automatically determined points within tgt pattern

# Holey support settings
tgtPattern      = False     # use same tgt pattern on different stage positions (useful for collection on holey support film)
alignToP        = False     # use generic View image in buffer P to align to
refineVec       = False     # refine tgt pattern for local stage position by aligning furthest targets along each axis to to buffer P
refineGeo       = False     # uses on-the-fly CtfFind results of first image to refine geometry before tilting (only use when CTF fits on your sample seem reliable)

# Session settings
beamTiltComp    = True      # use beam tilt compensation (uses coma vs image shift calibrations)
addAF           = False     # does autofocus at the start of every tilt group, increases exposure on tracking TS drastically
previewAli      = False      # adds initial dose, but makes sure start tilt image is on target (uses view image and aligns to buffer P if alignToP == True)
viewAli         = True     # adds an alignment step with a View image if it was saved during the target selection (only if previewAli is activated)

# Output settings
sortByTilt      = True      # sorts tilt series by tilt angle after acquisition is completed (takes additional time), requires mrcfile module
binFinalStack   = 1         # bin factor for final saved stack after acquisition (unbinned stack will be deleted to save storage space)
delFinalStack   = False     # delete final tilt series stacks (only keeps frames for reconstruction to save storage space) 
doCtfFind       = False     # set to False to skip CTFfind estimation (only necessary if it causes crashes => if it does crash, SerialEM will output some troubleshoot data that you should send to David!) 
doCtfPlotter    = True      # runs ctfplotter instead of CTFfind, needs standalone version of 3dmod on PATH
extendedMdoc    = True      # saves additional info to .mdoc file

# Hardware settings
slowTilt        = False     # do backlash step for all tilt angles, on bad stages large tilt steps are less accurate
fixedStageTilt  = False     # keep stage at fixedStageTiltAngle while running the full scheduled tilt series
fixedStageTiltAngle = 0.0   # physical stage angle [degrees] when fixedStageTilt is True
taOffsetPos     = -1.71         # additional tilt axis offset values [microns] applied to calculations for positive and...
taOffsetNeg     = -1.71         # ...negative branch of the tilt series (possibly useful for side-entry holder systems)
checkDewar      = True      # check if dewars are refilling before every acquisition
cryoARM         = False     # if you use a JEOL cryoARM TEM, this will keep the dewar refilling in sync
coldFEG         = True     # if you use a cold FEG, this will flash the gun whenever the dewars are being refilled
flashInterval   = 0        # time in hours between cold FEG flashes, -1: flash only during dewar refill (interval is ignored on Krios, uses FlashingAdvised function instead)
flashErrorWaitTime = 180    # wait time if error occurs in flashing
slitInterval    = 0         # time in minutes between centering the energy filter slit using RefineZLP, ONLY works with tgtPattern (needs pattern vectors to find good position for alignment)

# Advanced settings
fitLimit        = 30        # refineGeo: minimum resolution [Angstroms] needed for CTF fit to be considered for refineGeo
parabolTh       = 9         # refineGeo: minimum number of passable CtfFind values to fit paraboloid instead of plane 
imageShiftLimit = 20        # maximum image shift [microns] SerialEM is allowed to apply (this is a SerialEM property entry, default is 15 microns)
dataPoints      = 4         # number of recent specimen shift data points used for estimation of eucentric offset (default: 4)
alignLimit      = 0.5       # maximum shift [microns] allowed for record tracking between tilts, should reduce loss of target in case of low contrast (not applied for tracking TS); also the threshold to take a second tracking image when using trackTwice
minCounts       = 0         # minimum mean counts per second of record image (if set > 0, tilt series branch will be aborted if mean counts are not sufficient)
ignoreNegStart  = True      # ignore first shift on 2nd branch, which is usually very large on bad stages
realignToItem   = False     # Use SerialEM's RealignToItem routine instead of simple image realignment (was default in PACEtomo <=v1.9.1)
refFromPreview  = False     # Makes temporary reference from Preview image collected during previewAli for use with first Record image
noZeroRecAli    = False     # Skip alignment of first tilt image to reference 
autoStartTilt   = False     # Uses measured pretilt to set compensating startTilt      
tiltTargets     = 0         # Stage tilt at which targets were selected (if not 0, it will be automatically used as startTilt!)

# Target montage settings
tgtMontage      = False     # collect montage for each target using the shorter camera dimension (e.g. for square aperture montage tomography)
tgtMntSize      = 1         # size of montage pattern (1: 3x3, 2: 5x5, 3: 7x7, ...)
tgtMntOverlap   = 0.05      # montage tile overlap as fraction of shorter camera dimension
tgtMntXOffset   = 0         # max offset [microns] applied along tilt axis throughout tilt series (+tgtMntXOffset is reached at maxTilt, -tgtMntXOffset at minTilt)
tgtMntFocusCor  = False     # do focus compensation for tiles of montage
tgtTrackMnt     = False     # set to True if you also want the tracking target to be a montage

debug           = False     # Enables additional output for a few processes (e.g. cross-correlation for all image alignments)
breakpoints     = False     # Waits at every debug output for user to press B key.

# Defocus measure / autofocus (replaces sem.G / sem.G(-1))
# ctf | beam_tilt (physics+calibration) | beam_tilt_sem (sem.G(-1)+calibration)
defocusMethod = "beam_tilt_sem"
# False on scopes without XLensDeflector: skip all XLens Report/Set/Restore.
# Requires doRonchigram = False.
hasXLens = True
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
ctfDefocusLo = -12.0            # CtfFind search range low [microns]
ctfDefocusHi = -0.2             # CtfFind search range high [microns]
ctf_resolution_max_A = 20.0     # retry CtfFind if resolution [A] is above this
ctf_max_attempts = 3
ctf_retry_delay_s = 5
tilt_angle_mrad = 10.0            # match calibrate_beam_tilt_scaling.py
beam_tilt_correction = None       # SetBeamTilt scale; same value used in defocus beta
defocus_tilt_correction = beam_tilt_correction
beam_tilt_xtilt_x = 0.0           # X-tilt for beam-tilt defocus (ignored when hasXLens is False)
beam_tilt_xtilt_y = 0.0
spherical_aberration_mm = 2.7     # Cs for defocus = -disp/(2*beta) - Cs*beta^2 [mm]
autofocus_cycles = 2
measure_cycles = 1
autofocus_tolerance_um = 0.05

########## Ronchigram / laser alignment ##########
# Trial LD area must match Record position; only exposure should differ.
# Overridable from target file via _bset (e.g. _bset doRonchigram true).
# Requires hasXLens = True.
doRonchigram       = True
ronchiBaseSuffix   = "_ronchi"         # appended to active frame base name for Trial saves only, then restored
ronchiC3Offset     = -20          # added to ReportImageDistanceOffset before Trial shot
ronchiDelay        = 1.0          # seconds after C3 offset change
ronchiBinning      = 32
ronchiPixelSize    = 0.98e-4 * 2 # um (unbinned; multiplied by binning in analysis)
ronchiTargetPhaseA = -1.93941993           # vertical laser (rad)
ronchiTargetPhaseB = 1.67658165        # horizontal laser (rad)
ronchiCorrectKs    = [[9.303, -0.662] ,  [0.856 ,8.680]]
ronchiPeakRadius   = 100
ronchiMontage      = True         # also run before montage tile Record shots
ronchiCorrMatrix   = [[0.212, 1.28], [1.22, -0.243]]  # phase-to-deflector coupling, scaled by 1e-5
ronchiCorrectC3    = True         # apply C3 correction from mean ks error (diagonal fringe spacing)
ronchiC3CorrectionFactor = 20 / 9.1  # um offset per um^-1 mean ks error
ronchiMinErrForC3Correction   = 0.3          # apply C3 on 1st Trial only if |c3 correction| exceeds this (um)
ronchiMinErrForC3CorrectionRedo = 0.5        # apply C3 on 2nd Trial only if |c3 correction| exceeds this (um)
redo_ronchi_after_C3 = True       # up to 3 Trials: 1st C3, 2nd optional C3 + 3rd phase-only if 2nd C3 applied
ronchiPerPositionC3 = True        # remember ImageDistanceOffset per target; False = global C3 for all
ronchiXLensTolerance = 0.000125     # reset XLensDeflector(2) to start if |x-x0| or |y-y0| exceeds this
ronchiStartXLensX = None          # set from ReportXLensDeflector(2) at startup when doRonchigram
ronchiStartXLensY = None
ronchiStartC3Offset = None      # set from ReportImageDistanceOffset at startup when doRonchigram
########## END Ronchigram settings ##########

######### LAFIS: lpp afis correction #########
# calibration matrix applied when beamTiltComp == True on xlpp
# Requires hasXLens = True and beamTiltComp = True to be meaningful.

xt_is_matrix = [[0.000324, -0.000347],[0.001100, 0.00028125]]  #26jul23
df_is_matrix = [[0.041381,0.012342], [0.041381,0.012342]]

lafisZeroImageShiftDefocus = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisZeroImageShiftXLens = None            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisIsDone = False            # set from saveZeroImageShiftDefocusXLens before doLafis
lafisXtCorrectionX = 0.0       # set from doLafis as the correction made on XLens
lafisXtCorrectionY = 0.0       # set from doLafis as the correction made on XLens
########## END Lafis settings ##########

########## END SETTINGS ########## 

default_startTilt = startTilt
default_minTilt = minTilt
default_maxTilt = maxTilt
default_pretilt = pretilt
default_rotation = rotation

versionPACE = "1.9.2c"
import sys
sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
import serialem as sem
import os
import copy
import time
import struct
from datetime import datetime
import glob
from functools import wraps
import numpy as np
from scipy import optimize
import PACEtomo_beamTiltDefocus as btdef

if sortByTilt: 
    import subprocess
    import mrcfile
    from pathlib import Path

btdef.configure(sem_module=sem, has_x_lens=hasXLens)

# Per-run session state (initialized in run_one_nav_item; module defaults for helpers/linter)
geo = [[], [], []]
posResumed = -1
resumePN = 0
resumePlus = 0
resumeMinus = 0
resumePercent = 0
maxProgress = 1
startTime = 0
tiltStepCounter = 0
schedule_start = 0
is2ssMatrix = ss2isMatrix = s2ssMatrix = c2ssMatrix = ss2cMatrix = None
camX = camY = 0
origMag = 0

versionCheck = sem.IsVersionAtLeast("40200", "20240814")
if not versionCheck and sem.IsVariableDefined("warningVersion") == 0:
    runScript = sem.YesNoBox("\n".join(["WARNING: You are using a version of SerialEM that does not support all PACEtomo features. It is recommended to update to the latest SerialEM beta version!", "", "Do you want to run PACEtomo regardless?"]))
    if not runScript:
        sem.Exit()
    else:
        sem.SetPersistentVar("warningVersion", "")

########### FUNCTIONS ###########

def checkFilling():
    global dewarFillTime
    filling = sem.AreDewarsFilling()
    timerStart = 0
    if filling >= 1:
        timerStart = time.time()
        log(datetime.now().strftime("%d.%m.%Y %H:%M:%S") + ": Dewars are filling...", color=4)
        if cryoARM:                                                                             # make sure both tanks are being filled on cryoARM
            sem.LongOperation("RS", "0", "RT", "0")
        if coldFEG:                                                                             # flash gun while dewars refill
            checkColdFEG()
    while filling >= 1:
        log("Dewars are still filling...")
        sem.Delay(60, "s")
        filling = sem.AreDewarsFilling()
    if timerStart > 0:
        log(f"Dewars finished filling after {(time.time() - timerStart) / 60} minutes.")
        dewarFillTime = time.time() - timerStart

def checkColdFEG():
    if not cryoARM:                                                                             # Routine for Krios CFEG with Advanced scripting >4
        flashLow = sem.IsFEGFlashingAdvised(0)
        flashHigh = sem.IsFEGFlashingAdvised(1)
        if flashHigh == 1:
            sem.NextFEGFlashHighTemp(1)
        else:
            flashLow = sem.IsFEGFlashingAdvised(0)
            sem.NextFEGFlashHighTemp(0)
        if flashLow == 1 or flashHigh ==1:
            try:
                log("Flashing cold feg....")
                sem.LongOperation("FF", "0")
            except Exception as e:
                time.sleep(flashErrorWaitTime)
            log("Flashing cold feg done")
    else:
            sem.LongOperation("FF", str(flashInterval))

def checkSlit(vec, size, tilt, pn):                                                             # check ZLP in hole outside of pattern along tilt axis
    global lastSlitCheck
    log("Refining ZLP...", style=1)
    sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
    shift = vec * (size + 1)
    shift[1] *= np.cos(np.radians(tilt))
    sem.ImageShiftByMicrons(*shift)
    sem.RefineZLP()
    sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
    lastSlitCheck = sem.ReportClock()

def checkValves():
    if not int(sem.ReportColumnOrGunValve()):
        sem.SetColumnOrGunValve(1)


def ensure_fixed_stage_tilt(when="", quiet=False):
    """Move stage to fixedStageTiltAngle when fixedStageTilt is enabled."""
    if not fixedStageTilt:
        return
    target = float(fixedStageTiltAngle)
    current = float(sem.ReportTiltAngle())
    if abs(current - target) > 0.05:
        sem.TiltTo(target)
        if not quiet:
            suffix = f" {when}" if when else ""
            log(f"NOTE: fixedStageTilt: stage set to {target:.1f} deg{suffix}.")


def serialEM_measure_defocus():
    """SerialEM sem.G(-1) defocus with optional scaling calibration."""
    xt_x = beam_tilt_xtilt_x if hasXLens else None
    xt_y = beam_tilt_xtilt_y if hasXLens else None
    return btdef.measure_serialEM_defocus(
        xtilt_x=xt_x,
        xtilt_y=xt_y,
    )


def beam_tilt_measure_defocus():
    """Beam-tilt defocus via shared calibrated measurement."""
    xt_x = beam_tilt_xtilt_x if hasXLens else None
    xt_y = beam_tilt_xtilt_y if hasXLens else None
    return btdef.measure_defocus(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        defocus_tilt_correction=defocus_tilt_correction,
        xtilt_x=xt_x,
        xtilt_y=xt_y,
        cs_mm=spherical_aberration_mm,
    )


def ctf_measure_defocus():
    """CTF defocus: set X-tilt, Focus, CtfFind (with retries), restore X-tilt."""
    xtX = xtY = None
    if hasXLens:
        xtX, xtY = sem.ReportXLensDeflector(2)
    try:
        if hasXLens:
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
                    log(f"ERROR: CtfFind failed on attempt {attempt}/{ctf_max_attempts}.")
                    if attempt < ctf_max_attempts:
                        continue
                    return np.nan
                resolution = float(cfind[-1])
                if resolution <= ctf_resolution_max_A:
                    break
                log(
                    f"WARNING: CtfFind resolution {resolution:.2f} A > {ctf_resolution_max_A} A "
                    f"(attempt {attempt}/{ctf_max_attempts})"
                )
            else:
                log(
                    f"WARNING: CtfFind resolution still > {ctf_resolution_max_A} A after "
                    f"{ctf_max_attempts} attempts; using last result."
                )
        finally:
            sem.NoMessageBoxOnError(0)
        if len(cfind) == 0:
            return np.nan
        defocus = float(cfind[0])
        log(f"CtfFind: {defocus:.4f} microns ({float(cfind[-1]):.2f} A)")
        return defocus
    finally:
        if hasXLens and xtX is not None:
            sem.SetXLensDeflector(2, xtX, xtY)


def measure_defocus():
    """sem.G(-1): measure defocus only, no focus change."""
    if defocusMethod not in ("ctf", "beam_tilt", "beam_tilt_sem"):
        sem.OKBox(
            f"ERROR: Unknown defocusMethod '{defocusMethod}'. "
            "Use 'ctf', 'beam_tilt', or 'beam_tilt_sem'."
        )
        sem.Exit()
    if defocusMethod == "beam_tilt":
        defocus = np.nan
        speed_x = speed_y = 0.0
        for _ in range(measure_cycles):
            defocus, speed_x, speed_y = beam_tilt_measure_defocus()
        return float(defocus), np.array([speed_x, speed_y])
    if defocusMethod == "beam_tilt_sem":
        defocus = np.nan
        for _ in range(measure_cycles):
            defocus, speed_x, speed_y = serialEM_measure_defocus()
        return float(defocus), np.array([speed_x, speed_y])
    defocus = ctf_measure_defocus()
    return float(defocus), np.array([0.0, 0.0])


def autofocus_apply(target):
    """sem.G: measure defocus and correct to target."""
    defocus = np.nan
    for cycle in range(1, autofocus_cycles + 1):
        if defocusMethod == "beam_tilt":
            defocus, speed_x, speed_y = beam_tilt_measure_defocus()
        elif defocusMethod == "beam_tilt_sem":
            defocus, speed_x, speed_y = serialEM_measure_defocus()
        else:
            defocus = ctf_measure_defocus()
            speed_x = speed_y = 0.0
        if not np.isfinite(defocus):
            log(f"WARNING: Autofocus measurement failed on cycle {cycle}/{autofocus_cycles}.")
            return defocus
        error = target - defocus
        log(
            f"Autofocus {cycle}/{autofocus_cycles}: measured={defocus:.4f} um, "
            f"target={target:.3f} um, error={error:.3f} um"
        )
        if abs(error) <= autofocus_tolerance_um:
            return defocus
        sem.ChangeFocus(error)
    return defocus

##############################################################################
# Ronchigram - image analysis (numpy / FFT only; no SerialEM calls)
#
# Public API: analyze_ronchigram()
# Internal helpers: _ronchi_bin_image, _ronchi_find_fourier_centered, ...
##############################################################################


def _ronchi_bin_image(image, binning=32):
    bins = [image[i:(image.shape[0] // binning * binning):binning, j:(image.shape[1] // binning * binning):binning]
            for i in range(binning) for j in range(binning)]
    return np.sum(bins, axis=0)


def _ronchi_find_fourier_centered(image, padded_size=4096):
    fourier = np.fft.fftshift(np.fft.fft2(image - np.mean(image), s=(padded_size, padded_size)))
    center = np.array(image.shape) / 2
    grid_y, grid_x = np.indices((padded_size, padded_size), dtype=float)
    shifted_y = grid_y - padded_size / 2
    shifted_x = grid_x - padded_size / 2
    correction_phase_x = np.exp(2j * np.pi * shifted_x / padded_size * center[1])
    correction_phase_y = np.exp(2j * np.pi * shifted_y / padded_size * center[0])
    return fourier * correction_phase_x * correction_phase_y


def _ronchi_report_angles(ks, start_angle=-135):
    return np.mod(np.arctan2(-ks[:, 0], ks[:, 1]) * 180 / np.pi - start_angle, 360) + start_angle


def _ronchi_find_peaks(fourier, radius=100, npeaks=4):
    fourier_abs = np.abs(fourier).copy()
    peak_locations = []
    phases = []
    size = np.shape(fourier_abs)[0]
    peak_coords = [size // 2, size // 2]
    fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0
    for _ in range(npeaks):
        max_idx = np.argmax(fourier_abs)
        peak_coords = np.unravel_index(max_idx, fourier_abs.shape)
        peak_locations.append(peak_coords)
        fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                    (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0
        phases.append(np.angle(fourier[peak_coords]))
    return np.array(peak_locations) - size / 2, np.array(phases)


def _ronchi_find_ks_phases(corrected_fourier, pixel_size_um, npeaks=2, radius=100, binning=1, fourier_size=None):
    if fourier_size is None:
        fourier_size = corrected_fourier.shape[0]
    peaks, phases = _ronchi_find_peaks(corrected_fourier, radius=radius, npeaks=npeaks * 2)
    ordering = np.argsort(_ronchi_report_angles(peaks, start_angle=-135))
    peaks = peaks[ordering][:npeaks]
    phases = phases[ordering][:npeaks]
    ks = peaks / fourier_size * 1 / (pixel_size_um * binning)
    return ks, phases


def analyze_ronchigram(image, pixel_size_um, binning, target_phase_a, target_phase_b,
                       correct_ks, peak_radius=100, corr_matrix=None, corr_scale=1e-5,
                       c3_correction_factor=20 / 6.85):
    """FFT peak phases -> laser deflector and C3 corrections (dict with correction_x/y, phases, ks, ...)."""
    if corr_matrix is None:
        corr_matrix = [[0.212, 1.28], [1.22, -0.243]]
    correct_ks = np.asarray(correct_ks, dtype=float)
    binned = _ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = _ronchi_find_fourier_centered(binned)
    ks, phases = _ronchi_find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius, binning=1,
                                       fourier_size=image_fft.shape[0])
    ks_error = ks - correct_ks
    ks_total_err = float(np.linalg.norm(ks_error))
    ks_avg_err = (ks_error[0, 0] + ks_error[1, 1]) / 2.0
    c3_correction = ks_avg_err * c3_correction_factor
    phase_err_a = np.mod(phases[0] - target_phase_a + np.pi, 2 * np.pi) - np.pi
    phase_err_b = np.mod(phases[1] - target_phase_b + np.pi, 2 * np.pi) - np.pi
    corr = np.asarray(corr_matrix, dtype=float) * corr_scale
    correction_x = phase_err_a * corr[0, 0] + phase_err_b * corr[0, 1]
    correction_y = phase_err_a * corr[1, 0] + phase_err_b * corr[1, 1]
    return {
        "ks": ks, "phases": phases, "ks_error": ks_error, "ks_total_err": ks_total_err,
        "ks_avg_err": ks_avg_err,
        "c3_correction": c3_correction,
        "phase_err_a": phase_err_a, "phase_err_b": phase_err_b,
        "correction_x": correction_x, "correction_y": correction_y,
    }


##############################################################################
# Ronchigram - microscope aligner (SerialEM: Trial acquire, apply, Record)
#
# Public API: checkRonchigramSetup(), doRonchigramCorrection(), recordWithRonchi()
##############################################################################


def applyRonchigramXtiltCorrection(correction_x, correction_y, lens_index=2):
    """Apply calculated shifts to the X lens deflector."""
    if not hasXLens:
        return
    xtX, xtY = sem.ReportXLensDeflector(lens_index)
    sem.SetXLensDeflector(lens_index, xtX + correction_x, xtY + correction_y)


def _reset_ronchi_xlens_if_out_of_tolerance(lens_index=2):
    """Reset XLensDeflector to session start if phase corrections have walked it too far."""
    global ronchiStartXLensX, ronchiStartXLensY
    if not hasXLens or ronchiStartXLensX is None or ronchiStartXLensY is None:
        return
    xtX, xtY = sem.ReportXLensDeflector(lens_index)
    xtX, xtY = float(xtX), float(xtY)
    xtX_delta = xtX - lafisXtCorrectionX - ronchiStartXLensX
    xtY_delta = xtY - lafisXtCorrectionY - ronchiStartXLensY
    if (abs(xtX_delta) > ronchiXLensTolerance
            or abs(xtY_delta) > ronchiXLensTolerance):
        log(
            f"WARNING: Ronchigram X lens deflector ({xtX:.6f}, {xtY:.6f}) beyond tolerance "
            f"{ronchiXLensTolerance} from start ({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}); resetting."
        )
        sem.SetXLensDeflector(lens_index, ronchiStartXLensX+lafisXtCorrectionX, ronchiStartXLensY+lafisXtCorrectionY)


def applyRonchigramC3Correction(c3_correction, baseline_offset):
    """Apply C3 correction to ImageDistanceOffset relative to pre-ronchigram offset."""
    new_offset = baseline_offset + c3_correction
    sem.SetImageDistanceOffset(new_offset)
    return new_offset


def _report_frame_basename():
    """ReportFrameBaseName -> (use_in_frame, name, use_in_folder)."""
    r = sem.ReportFrameBaseName()
    use_in_frame = int(r[0])
    name = r[1] if len(r) > 1 else ""
    use_in_folder = int(r[2]) if len(r) > 2 else 0
    if isinstance(name, str) and name.lower() == "none":
        name = ""
    return use_in_frame, name, use_in_folder


def _ronchi_trial_basename(use_in_frame, current_name, use_in_folder):
    """Frame base name for Trial ronchigram image (active SetFrameBaseName + suffix)."""
    root = (current_name or "").strip()
    return root + ronchiBaseSuffix, root


def _set_ronchi_trial_frame_basename():
    """SetFrameBaseName with _ronchi suffix for Trial; return saved state for restore."""
    use_in_frame, name, use_in_folder = _report_frame_basename()
    ronchi_name, root = _ronchi_trial_basename(use_in_frame, name, use_in_folder)
    if not root:
        log("WARNING: Ronchigram Trial: no frame base name; SetFrameBaseName before acquire.")
    else:
        sem.SetFrameBaseName(0, use_in_frame, use_in_folder, ronchi_name)
        log(f"Ronchigram Trial: SetFrameBaseName -> {ronchi_name}")
    return use_in_frame, name, use_in_folder


def _restore_frame_basename(saved):
    use_in_frame, name, use_in_folder = saved
    sem.SetFrameBaseName(0, use_in_frame, use_in_folder, name)
    log(f"Ronchigram Trial: restored frame base name -> {name or '(none)'}")


def checkRonchigramSetup():
    """Startup checks when doRonchigram is enabled (called after parseTargets)."""
    global ronchiStartXLensX, ronchiStartXLensY, ronchiStartC3Offset, doRonchigram
    if not hasXLens and doRonchigram:
        sem.OKBox(
            "ERROR: doRonchigram requires XLensDeflector. "
            "Set hasXLens = True, or set doRonchigram = False."
        )
        sem.Exit()
    if not doRonchigram:
        if not hasXLens:
            log("NOTE: hasXLens=False; all XLens Report/Set/Restore calls are skipped.")
        return
    if ronchiStartXLensX is None:
        ronchiStartXLensX, ronchiStartXLensY = [float(v) for v in sem.ReportXLensDeflector(2)[:2]]
    if ronchiStartC3Offset is None:
        ronchiStartC3Offset = float(sem.ReportImageDistanceOffset())
    log(
        f"NOTE: Ronchigram X lens deflector start ({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}), "
        f"reset tolerance {ronchiXLensTolerance}"
    )
    log(f"NOTE: Ronchigram C3 (ImageDistanceOffset) session start {ronchiStartC3Offset:.2f} um")
    trial_exp, *_ = sem.ReportExposure("T")
    record_exp, *_ = sem.ReportExposure("R")
    if trial_exp <= 0 and sem.IsVariableDefined("warningRonchiTrial") == 0:
        sem.Pause("WARNING: Trial exposure is zero or not set. Configure Trial with a very short exposure at the same position as Record.")
        sem.SetPersistentVar("warningRonchiTrial", "")
    if trial_exp >= record_exp * 0.5:
        log("WARNING: Trial exposure should be much shorter than Record for negligible ronchigram dose.")
    try:
        t_shift = sem.ReportLDAreaShift("T")
        if len(t_shift) >= 2 and np.linalg.norm(np.array(t_shift[:2], dtype=float)) > 0.01:
            log("WARNING: Trial area position differs from Record. Set Trial LD offsets to match Record.")
    except (AttributeError, TypeError, ValueError):
        pass
    log("NOTE: Ronchigram uses Trial at Record beam position. Set Trial LD offsets identical to Record; only exposure should differ.")
    log(f"NOTE: Ronchigram Trial temporarily appends '{ronchiBaseSuffix}' to the active frame base name.")
    log("NOTE: Ronchigram runs before startup preview alignment images when previewAli or viewAli is enabled.")
    if ronchiCorrectC3:
        log(
            f"NOTE: Ronchigram C3 correction enabled (factor {ronchiC3CorrectionFactor:.4f}, "
            f"only if |C3 correction| > {ronchiMinErrForC3Correction} um)."
        )
    if redo_ronchi_after_C3:
        log(
            "NOTE: Ronchigram may use up to 3 Trials before Record: "
            "1st C3 (if above threshold), 2nd phase or C3, 3rd phase-only if 2nd C3 applied."
        )
        log(
            f"NOTE: 2nd-Trial C3 threshold |correction| > {ronchiMinErrForC3CorrectionRedo} um "
            f"(1st-Trial threshold {ronchiMinErrForC3Correction} um)."
        )
    if ronchiPerPositionC3:
        log("NOTE: Ronchigram C3 offset is stored and restored per target across tilts.")


def _apply_stored_c3_offset(pos, pn):
    """Restore per-target C3 offset, or initialize from current scope value."""
    global position
    stored = position[pos][pn].get("c3_offset")
    if stored is not None:
        sem.SetImageDistanceOffset(stored)
        log(f"Ronchigram: target {pos + 1} using stored C3 offset {stored:.2f} um")
        return float(stored)
    current = float(sem.ReportImageDistanceOffset())
    for b in range(len(position[pos])):
        position[pos][b]["c3_offset"] = current
    log(f"Ronchigram: target {pos + 1} initializing C3 offset {current:.2f} um")
    return current


def _save_c3_offset_for_target(pos, pn):
    """Save current ImageDistanceOffset for this target (all branches)."""
    global position
    current = float(sem.ReportImageDistanceOffset())
    for b in range(len(position[pos])):
        position[pos][b]["c3_offset"] = current
    log(f"Ronchigram: target {pos + 1} saved C3 offset {current:.2f} um")


def _log_ronchi_ks(result, pass_label=""):
    prefix = f"Ronchigram{pass_label}"
    log(
        f"{prefix} ks (1/um): {np.array2string(result['ks'], precision=4)} | "
        f"ks error: {np.array2string(result['ks_error'], precision=4)}"
    )
    log(
        f"{prefix} ||ks error||: {result['ks_total_err']:.4f} (1/um) | "
        f"mean diagonal ks error: {result['ks_avg_err']:.4f} (1/um) | "
        f"recommended C3 correction: {result['c3_correction']:.2f} um"
    )


def _log_ronchi_phases(result, pass_label=""):
    phases = result["phases"]
    prefix = f"Ronchigram{pass_label}"
    log(
        f"{prefix} phases (rad): measured vertical={phases[0]:.3f} horizontal={phases[1]:.3f} | "
        f"targets vertical={ronchiTargetPhaseA:.3f} horizontal={ronchiTargetPhaseB:.3f}"
    )
    log(
        f"{prefix} phase error (rad): vertical={result['phase_err_a']:.3f} "
        f"horizontal={result['phase_err_b']:.3f} | "
        f"deflector dX={result['correction_x']:.3e} dY={result['correction_y']:.3e}"
    )


def _analyze_ronchi_image(image):
    return analyze_ronchigram(
        image, ronchiPixelSize, ronchiBinning, ronchiTargetPhaseA, ronchiTargetPhaseB,
        ronchiCorrectKs, peak_radius=ronchiPeakRadius, corr_matrix=ronchiCorrMatrix,
        c3_correction_factor=ronchiC3CorrectionFactor,
    )


def _acquire_ronchi_trial(trial_offset_baseline, pass_label=""):
    """Trial ronchigram at C3 imaging offset; restore baseline offset after shot."""
    _reset_ronchi_xlens_if_out_of_tolerance()
    is_x, is_y, *_ = sem.ReportImageShift()
    sem.GoToLowDoseArea("T")
    sem.SetImageShift(0, 0)
    sem.SetImageShift(is_x, is_y)
    sem.SetImageDistanceOffset(trial_offset_baseline + ronchiC3Offset)
    saved_basename = _set_ronchi_trial_frame_basename()
    try:
        sem.Delay(ronchiDelay, "s")
        add_lpp_meta_to_next_mdoc()
        sem.T()
    finally:
        sem.SetImageDistanceOffset(trial_offset_baseline)
        _restore_frame_basename(saved_basename)
    if pass_label:
        log(f"Ronchigram{pass_label}: Trial image acquired.")
    return np.asarray(sem.bufferImage("A"))


def _try_apply_ronchi_c3(result, c3_baseline_offset, pass_label="", min_err=None):
    """Apply C3 if enabled and correction magnitude exceeds minimum. Returns True if C3 was changed."""
    if min_err is None:
        min_err = ronchiMinErrForC3Correction
    prefix = f"Ronchigram{pass_label}"
    if not ronchiCorrectC3:
        log(f"{prefix} C3: correction {result['c3_correction']:.2f} um not applied (ronchiCorrectC3=False)")
        return False
    if abs(result["c3_correction"]) <= min_err:
        log(
            f"{prefix} C3: skipped (|correction| {abs(result['c3_correction']):.2f} um <= "
            f"minimum {min_err} um)"
        )
        return False
    new_offset = applyRonchigramC3Correction(result["c3_correction"], c3_baseline_offset)
    log(
        f"{prefix} C3: ImageDistanceOffset adjusted by {result['c3_correction']:.2f} um "
        f"(now {new_offset:.2f} um)"
    )
    return True


def _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=""):
    """Acquire Trial ronchigram and analyze. Returns analysis result dict."""
    image = _acquire_ronchi_trial(c3_baseline_offset, pass_label=pass_label)
    result = _analyze_ronchi_image(image)
    _log_ronchi_ks(result, pass_label=pass_label)
    return result


def _apply_ronchi_phase(result, pass_label=""):
    _log_ronchi_phases(result, pass_label=pass_label)
    applyRonchigramXtiltCorrection(result["correction_x"], result["correction_y"])

def calc_xt_is(xt0, is_delta):
   xt1 = [0.0,0.0]
   xt1[0] = xt0[0]+is_delta[0]*xt_is_matrix[0][0]+is_delta[1]*xt_is_matrix[1][0]
   xt1[1] = xt0[1]+is_delta[0]*xt_is_matrix[0][1]+is_delta[1]*xt_is_matrix[1][1]
   return xt1

def calc_df_is(df0, is_delta):
   df1 = 0.0
   df1 =df0+is_delta[0]*df_is_matrix[0][0]+is_delta[1]*df_is_matrix[1][1]
   return df1

def add_lpp_meta_to_next_mdoc():
    for k,v in (
            ('ImageDistanceOffset', sem.ReportImageDistanceOffset()),
        ):
        v_str = '%.12f' % (float(v))
        sem.AddToNextFrameStackMdoc(k, v_str)

    for k,v in (
            ('BeamTilt', sem.ReportBeamTilt()),
            ('ObjectiveStig',sem.ReportObjectiveStigmator()),
            ('XTilt', sem.ReportXLensDeflector(2)),):
        v_str = '%.6f        %.6f' % (float(v[0]), float(v[1]))
        sem.AddToNextFrameStackMdoc(k, v_str)

def saveZeroImageShiftDefocusXLens():
    global lafisZeroImageShiftDefocus
    global lafisZeroImageShiftXLens
    lafisZeroImageShiftDefocus = sem.ReportDefocus()
    if hasXLens:
        lafisZeroImageShiftXLens = sem.ReportXLensDeflector(2)
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
    if hasXLens:
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
    if hasXLens and lafisZeroImageShiftXLens is not None:
        xt_x, xt_y = lafisZeroImageShiftXLens
        sem.SetXLensDeflector(2, xt_x, xt_y)
        lafisXtCorrectionX = 0.0
        lafisXtCorrectionY = 0.0
    lafisIsDone = False
    log('WARNING: Lafis restored')

def doRonchigramCorrection(set_track_fn=None, pos=None, pn=None):
    """Trial shot + analyze_ronchigram + C3 and/or laser correction; return to Record area."""
    if not doRonchigram:
        return
    try:
        sem.UpdateLowDoseParams("T")
    except AttributeError:
        pass
    tilt = float(sem.ReportTiltAngle())
    trial_exp, *_ = sem.ReportExposure("T")
    log(
        f"Ronchigram: Trial acquire at tilt {tilt:.1f} deg | "
        f"C3 offset {ronchiC3Offset} um | binning {ronchiBinning} | Trial exposure {trial_exp:.4g} s"
    )
    if ronchiPerPositionC3 and pos is not None and pn is not None:
        c3_baseline_offset = _apply_stored_c3_offset(pos, pn)
    else:
        c3_baseline_offset = float(sem.ReportImageDistanceOffset())
    try:
        result = _ronchi_trial_and_analyze(c3_baseline_offset)
        c3_changed = _try_apply_ronchi_c3(result, c3_baseline_offset)

        if c3_changed and redo_ronchi_after_C3:
            log("Ronchigram: phase correction deferred until after 2nd Trial.")
            c3_baseline_offset = float(sem.ReportImageDistanceOffset())
            log("Ronchigram: 2nd Trial after 1st C3 change.")
            result = _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=" (2nd)")
            c3_changed_redo = _try_apply_ronchi_c3(
                result, c3_baseline_offset, pass_label=" (2nd)",
                min_err=ronchiMinErrForC3CorrectionRedo,
            )
            if c3_changed_redo:
                log("Ronchigram: phase correction deferred until after 3rd Trial.")
                c3_baseline_offset = float(sem.ReportImageDistanceOffset())
                log("Ronchigram: 3rd Trial (phase correction only).")
                result = _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=" (3rd)")
                _apply_ronchi_phase(result, pass_label=" (3rd)")
            else:
                _apply_ronchi_phase(result, pass_label=" (2nd)")
        else:
            _apply_ronchi_phase(result)

        if debug:
            log(f"DEBUG: Ronchigram pixel size {ronchiPixelSize} um, peak radius {ronchiPeakRadius} px")
        if ronchiPerPositionC3 and pos is not None and pn is not None:
            _save_c3_offset_for_target(pos, pn)
    except Exception as e:
        log(f"WARNING: Ronchigram analysis failed: {e}. Continuing without laser correction.")
        sem.SetImageDistanceOffset(ronchiStartC3Offset)
        if ronchiPerPositionC3 and pos is not None and pn is not None:
            _save_c3_offset_for_target(pos, pn)
    sem.GoToLowDoseArea("R")
    if set_track_fn is not None:
        set_track_fn()


def ronchi_before_preview_align(acquire_label="preview alignment", pos=None, pn=None):
    """Ronchigram laser correction at Record before Preview (L) alignment images."""
    if not doRonchigram:
        return
    log(f"NOTE: {acquire_label} - ronchigram before preview alignment image.")
    is_x, is_y, *_ = sem.ReportImageShift()
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    sem.SetImageShift(is_x, is_y)
    doRonchigramCorrection(set_track_fn=None, pos=pos, pn=pn)
    if beamTiltComp:
        log("Ronchigram: We should not resotre the value here?.")

def recordWithRonchi(set_track_fn=None, run_ronchi=True, acquire_label="Record", pos=None, pn=None):
    """Optional ronchigram correction, then standard Record acquire (sem.R / sem.S)."""
    if ronchiPerPositionC3 and pos is not None and pn is not None and not (run_ronchi and doRonchigram):
        _apply_stored_c3_offset(pos, pn)
    if run_ronchi and doRonchigram:
        log(f"NOTE: {acquire_label} - ronchigram Trial then Record stack frame.")
        doRonchigramCorrection(set_track_fn=set_track_fn, pos=pos, pn=pn)
    elif doRonchigram and not run_ronchi:
        log(f"NOTE: {acquire_label} - Record only (ronchigram skipped for this shot).")
    else:
        log(f"NOTE: {acquire_label} - Record acquire.")
    sem.Delay(delayIS, "s")
    add_lpp_meta_to_next_mdoc()
    sem.R()
    sem.S()

########### PACEtomo functions (non-ronchigram) ###########

def retryOpen(max_attempts=5, delay=5):
    """Decorator to retry function on permission exception."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            while attempts < max_attempts:
                try:
                    return func(*args, **kwargs)
                except PermissionError as e:
                    attempts += 1
                    if attempts == max_attempts:
                        raise e
                    log(f"WARNING: File [{args[0]}] could not be opened. Trying again... [{attempts}]")
                    time.sleep(delay)
            return None
        return wrapper
    return decorator

def openOldFile(file_path, max_attempts=5, delay=5):
    # Allow for multiple attempts to open file to prevent crash when file is being copied during access
    attempts = 0
    while attempts < max_attempts:
        try:
            sem.NoMessageBoxOnError(1)
            sem.OpenOldFile(file_path)
            # Break on success
            break
        except Exception as e:
            attempts += 1
            if attempts == max_attempts:
                raise e
            log(f"WARNING: File [{file_path}] could not be opened. Trying again... [{attempts}]")
            time.sleep(delay)
        finally:
            sem.NoMessageBoxOnError(0)

@retryOpen()
def parseTargets(file_path):
    """Reads targets file."""

    with open(file_path) as f:                                                                         # open last tgts or tgts_run file
        targetFile = f.readlines()

    targets = []
    geoPoints = []
    savedRun = []
    branch = None
    resume = {"sec": 0, "pos": 0}
    for line in targetFile:
        col = line.strip(os.linesep).split(" ")
        if col[0] == "": continue
        if line.startswith("_set") and len(col) == 4:
            if col[1] in globals():
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {col[3]}")
                globals()[col[1]] = float(col[3])
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")
        elif line.startswith("_bset") and len(col) == 4:
            if col[1] in globals():
                val = True if col[3].lower() in ["true", "yes", "y", "on"] else False
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {val}")
                globals()[col[1]] = val
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")                
        elif line.startswith("_spos"):
            resume["sec"] = int(col[2].split(",")[0])
            resume["pos"] = int(col[2].split(",")[1])
        elif line.startswith("_tgt"):
            targets.append({})
            branch = None
        elif line.startswith("_pbr"):
            savedRun.append([{},{}])
            branch = 0
        elif line.startswith("_nbr"):
            branch = 1
        elif line.startswith("_geo"):
            geoPoints.append({})
            branch = "geo"
        else:
            if branch is None:
                targets[-1][col[0]] = col[2]
            elif branch == "geo":
                geoPoints[-1][col[0]] = float(col[2])
            else:
                savedRun[-1][branch][col[0]] = col[2]
    if savedRun == []: savedRun = False
    return targets, savedRun, resume, geoPoints

@retryOpen()
def updateTargets(fileName, targets, position=[], sec=0, pos=0):
    output = ""
    if sec > 0 or pos > 0:
        output += "_set startTilt = " + str(startTilt) + "\n"
        output += "_set minTilt = " + str(minTilt) + "\n"
        output += "_set maxTilt = " + str(maxTilt) + "\n"
        output += "_set step = " + str(step) + "\n"
        output += "_set pretilt = " + str(pretilt) + "\n"
        output += "_set rotation = " + str(rotation) + "\n"
        output += "_spos = " + str(sec) + "," + str(pos) + "\n" * 2
    for pos in range(len(targets)):
        output += "_tgt = " + str(pos + 1).zfill(3) + "\n"
        for key in targets[pos].keys():
            output += key + " = " + targets[pos][key] + "\n"
        if position != []:
            output += "_pbr" + "\n"
            for key in position[pos][1].keys():
                output += key + " = " + str(position[pos][1][key]).strip("[]").replace("np.float64(", "").replace(")", "").replace(" ","") + "\n"   # temp fix for numpy>=2.0 changing str output for list of numbers to include np.float64()
            output += "_nbr" + "\n"
            for key in position[pos][2].keys():
                output += key + " = " + str(position[pos][2][key]).strip("[]").replace("np.float64(", "").replace(")", "").replace(" ","") + "\n"   # temp fix for numpy>=2.0 changing str output for list of numbers to include np.float64()     
        output += "\n"
    with open(fileName, "w") as f:
        f.write(output)

def geoPlane(x, a, b):
    return a * x[0] + b * x[1]

def geoPara(x, a, b, c, d, e):
    return a * x[0] + b * x[1] + c * (x[0]**2) + d * (x[1]**2) + e * x[0] * x[1]

def parseMdoc(mdocFile):
    with open(mdocFile, "r") as f:
        content = f.readlines()
    header = []
    items = []
    newItem = {}
    index = 0

    for line in content:
        col = line.strip().split()
        if len(col) == 0: 
            if "ZValue" in newItem.keys():
                items.append(newItem)
                newItem = {}
                continue
            else:
                continue
        if line.startswith("[ZValue"):
            index += 1
            newItem = {"index": index, "ZValue": col[2].strip("]")}
        elif "ZValue" in newItem.keys():
            newItem[col[0]] = [val for val in col[2:]]
        else:
            header.append(line.strip())
    if "ZValue" in newItem.keys(): # append last item
        items.append(newItem)
    return header, items

def writeMdoc(header, items, filename):
    content = ""
    for line in header:
        if line.startswith("[T"):
            content += os.linesep
        content += line + os.linesep
    content += os.linesep
    for i, item in enumerate(items):
        content += "[ZValue = " + str(i) + "]" + os.linesep
        item.pop("ZValue")
        item.pop("index")
        for key, attr in item.items():
            content += key + " = " + (" ".join(attr) if key != "DateTime" else "  ".join(attr)) # Date and time is separated by double space in SerialEM
            content += os.linesep
        if i < len(items) - 1: # Don't add additional linebreak at end of file
            content += os.linesep
    with open(filename, "w", newline="") as f:
        f.write(content)
    return

def getExtendedHeader(filename):
    with open(filename, "rb") as mrc_file:
        mrc = mrc_file.read(4096)

    # MRC header format SERI (SERI is not covered by mrcfile)
    format = ""
    for char in struct.iter_unpack("s", mrc[104:108]):
        format += char[0].decode("utf-8")
    log(f"DEBUG: MRC format: {format}")

    if format == "SERI":
        # Get number of sections
        section_number = struct.unpack("i", mrc[8: 12])[0]
        log(f"DEBUG: Number of sections: {section_number}")

        # Get bytes per section
        bytes_per_section = struct.unpack("h", mrc[128: 130])[0]
        log(f"DEBUG: Bytes per section: {bytes_per_section}")

        # Bitflags
        bitflags = struct.unpack("h", mrc[130: 132])[0]
        log(f"DEBUG: Bitflags: {bitflags}")
        """
        https://bio3d.colorado.edu/imod/doc/mrc_format.txt
        1 = Tilt angle in degrees * 100  (2 bytes)
        2 = X, Y, Z piece coordinates for montage (6 bytes)
        4 = X, Y stage position in microns * 25   (4 bytes)
        8 = Magnification / 100 (2 bytes)
        16 = Intensity * 25000  (2 bytes)
        32 = Exposure dose in e-/A2, a float in 4 bytes
        128, 512: Reserved for 4-byte items
        64, 256, 1024: Reserved for 2-byte items
        If the number of bytes implied by these flags does
        not add up to the value in nint, then nint and nreal
        are interpreted as ints and reals per section
        """

        section_data = []
        for i in range(1024, 1024 + bytes_per_section * section_number, bytes_per_section):     # extended header starts at byte 1024
            section_data.append(mrc[i:i + bytes_per_section])

        return section_data
    
    return None

def writeExtendedHeader(filename, section_data):
    with open(filename, "r+b") as mrc_file:
        mrc_file.seek(1024)
        mrc_file.write(b"".join(section_data))

@retryOpen()
def sortTS(ts_name):
    # Check for SPACEtomo installation to run sorting in background
    try:
        import SPACEtomo
        import inspect
        cli = Path(inspect.getfile(SPACEtomo)).parent / "CLI.py"
        SPACEtomo_version = SPACEtomo.__version__
    except ImportError:
        SPACEtomo_version = "0.0.0"
    
    from packaging.version import Version
    if Version(SPACEtomo_version) >= Version("1.3.1b19"):
        log(f"Sorting {ts_name} by tilt angle in background using SPACEtomo...")
        DETACHED_PROCESS = 0x00000008 # From here: https://stackoverflow.com/questions/89228/calling-an-external-command-in-python#2251026
        try:
            subprocess.Popen([sys.executable, cli, "sort", Path(curDir) / f"{ts_name}"], creationflags=DETACHED_PROCESS)
        except ValueError:      # Creationflags only supported on Windows
            subprocess.Popen([sys.executable, cli, "sort", Path(curDir) / f"{ts_name}"])
        return
    log(f"NOTE: Consider installing or updating SPACEtomo to sort tilt series in background!")

    log(f"Sorting {ts_name} by tilt angle...")
    if os.path.exists(os.path.join(curDir, ts_name + ".mdoc")):
        # Read mdoc file
        header, tilts = parseMdoc(os.path.join(curDir, ts_name + ".mdoc"))
        # Make list of tilt angles
        tiltAngles = [float(tilt["TiltAngle"][0]) for tilt in tilts]
        log(f"DEBUG: Tilts before sorting: {tiltAngles}")
        # Get extended header data
        ext_header_sections = getExtendedHeader(os.path.join(curDir, ts_name))
        log(f"DEBUG: MRC header sections: {len(ext_header_sections)}")
        # Open mrc file
        with mrcfile.open(os.path.join(curDir, ts_name), "r+") as mrc:
            stack = mrc.data
            log(f"DEBUG: MRC data dims: {np.array(stack).shape}")
            # Sort tilts and stack according to tilt angle
            zippedTilts = sorted(zip(tiltAngles, tilts, stack, ext_header_sections), key=lambda x: x[0])
            tiltAngles, tilts, stack, ext_header_sections = zip(*zippedTilts)
            log(f"DEBUG: Tilts after sorting: {tiltAngles}")
            stack = np.array(stack)
            # Save new stack in same file
            mrc.set_data(stack)
        # Update extended header
        writeExtendedHeader(os.path.join(curDir, ts_name), ext_header_sections)
        # Rename original mdoc
        os.rename(os.path.join(curDir, ts_name + ".mdoc"), os.path.join(curDir, os.path.splitext(ts_name)[0] + "_unsorted.mrc.mdoc"))
        # Save sorted mdoc
        writeMdoc(header, tilts, os.path.join(curDir, ts_name + ".mdoc"))
        log(f"NOTE: {ts_name} was sorted by tilt angle!")
    else:
        log(f"WARNING: {ts_name}.mdoc file could not be found! Tilt series stack was not sorted!")

def bin2d(img, factor):
    # Add third dimension if img is not stack
    if img.ndim == 2:
        img = np.expand_dims(img, 0)
    # Only bin in 2D and keep stack size
    factors = (1, factor, factor)
    # Calculate the new dimensions after cropping to allow even binning
    new_shape = tuple((dim // factor) * factor for dim, factor in zip(img.shape, factors))
    # Center crop the array to the new dimensions
    slices = tuple(slice((dim - new_dim) // 2, (dim + new_dim) // 2) for dim, new_dim in zip(img.shape, new_shape))
    cropped_img = img[slices]
    # Determine the new shape for reshaping
    reshaped_shape = np.array([(dim // factor, factor) for dim, factor in zip(cropped_img.shape, factors)]).reshape(-1)
    # Reshape the array
    reshaped_img = cropped_img.reshape(reshaped_shape)
    # Calculate the mean along the new axes
    for i in range(-1, -cropped_img.ndim-1, -1):
        reshaped_img = reshaped_img.mean(axis=i)
    # Remove added dimension
    if reshaped_img.shape[0] == 1:
        reshaped_img = reshaped_img[0, :, :]
    return reshaped_img

@retryOpen()
def binStack(ts_name, factor):
    factor = int(factor)
    log(f"Binning {ts_name} by {factor}...")
    if checkFrames(ts_name):
        # Check if tilt stack file exists
        if os.path.exists(os.path.join(curDir, ts_name)):
            # Open mrc file
            with mrcfile.open(os.path.join(curDir, ts_name), "r+") as mrc:
                stack = mrc.data
                voxel_size = mrc.voxel_size.x
                # Bin stack by factor
                stack = bin2d(stack, factor)
                # Save new stack in same file
                mrc.set_data(stack.astype(np.int16))
                # Update header pixel size
                mrc.voxel_size = voxel_size * factor
            log(f"NOTE: {ts_name} was binned by {factor}. Please use saved frames to regenerate the unbinned tilt series.")
        else:
            log(f"WARNING: {ts_name} file could not be found!")

def checkFrames(ts_name):
    # Check if frame saving is available by checking if warning was accepted by user at start of script
    if sem.IsVariableDefined("warningFramePath") == 0:
        # Make sure frames were saved
        frame_file, frame_dir, frame_name = sem.ReportLastFrameFile()
        if checkFrames and len(glob.glob(os.path.join(frame_dir, os.path.splitext(ts_name)[0] + "*"))) > 0:
            return True

    log(f"WARNING: Frames for {ts_name} could not be found. Keeping tilt stack unprocessed.")
    return False

def alignTo(buffer, debug=False):
    sem.AlignTo(buffer, 0, 0, 0, int(debug))
    if debug:
        try:
            sem.AddBufToStackWindow("A", 0, 0, 0, 0, "CC") #M #S [#B] [#O] [title]
        except AttributeError:
            # Show CC briefly, then switch back to aligned buffer for buffer shift
            sem.Delay(1, "s")
        sem.Copy("B", "A")
        sem.AlignTo(buffer)

def realignTo(nav_id=None, target=None):
    if target is not None and not realignToItem:
        # Move stage to target position
        sem.MoveStageTo(float(target["stageX"]), float(target["stageY"]))
        if "viewfile" in target.keys():
            sem.ReadOtherFile(0, "O", target["viewfile"]) # reads view file for first AlignTo instead
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            alignTo("O", debug)
            ASX, ASY = sem.ReportAlignShift()[4:6]
            log(f"Alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")    
        if "tgtfile" in target.keys():                
            sem.ReadOtherFile(0, "O", target["tgtfile"]) # reads tgt file for first AlignTo instead
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("R")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            if beamTiltComp:
                doLafis(is_x,is_y)
            ronchi_before_preview_align("initial realign preview (tgtfile)")
            sem.L()
            if beamTiltComp:
                restoreLafis()
            alignTo("O", debug)
            AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
            log(f"Alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
        elif "viewfile" in target.keys():
            # Use align between mags to align preview image to view image
            # If View image was already aligned, take new centered View image at startTilt and use as reference instead
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            sem.Copy("A", "O")
            # Check defocus offset
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("R") # Switch to R before applying defocus offset to not mess with potential mP/nP offsets between View and Rec
            sem.SetImageShift(0, 0)
            defocus_offset = max(-10, sem.ReportLDDefocusOffset("V"))
            if defocus_offset != 0:
                sem.ChangeFocus(defocus_offset) # Higher defocus for better correlation, but max at 10 to avoid major distortions
            sem.SetImageShift(is_x, is_y)
            if beamTiltComp:
                doLafis(is_x,is_y)
            ronchi_before_preview_align("initial realign preview (view to Record)")
            sem.L()
            if beamTiltComp:
                restoreLafis()
            sem.AlignBetweenMags("O", -1, -1, -1)
            AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
            if defocus_offset != 0:
                sem.ChangeFocus(-defocus_offset) # Reset focus
            log(f"Alignment (Pv2V) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
        else:
            log(f"WARNING: No target file or view file found for realignment!")
    elif nav_id is not None:
        sem.RealignToOtherItem(nav_id, 1)
    else:
        log(f"WARNING: No target provided for realignment!")

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

def breakpoint():
    """Breakpoint for debugging in SerialEM."""

    while not sem.KeyBreak():
        sem.Delay(0.1, "s")
    for i in range(5):
        if sem.KeyBreak("d"):
            dumpVars()
            break
        sem.Delay(0.1, "s")

def build_tilt_schedule(scheme, start_tilt, min_tilt, max_tilt, tilt_step, group_size):
    """Return ordered list of stage tilt angles [degrees] for the full series."""
    if scheme == "dose_symmetric":
        tilts = [start_tilt]
        plustilt = minustilt = start_tilt
        branchsteps = max(max_tilt - start_tilt, abs(min_tilt - start_tilt)) / group_size / tilt_step
        for _ in range(int(np.ceil(branchsteps))):
            for _ in range(group_size):
                plustilt += tilt_step
                if plustilt <= max_tilt + 1e-6:
                    tilts.append(plustilt)
            for _ in range(group_size):
                minustilt -= tilt_step
                if minustilt >= min_tilt - 1e-6:
                    tilts.append(minustilt)
        return tilts
    if scheme == "bidirectional":
        tilts = [start_tilt]
        angle = start_tilt + tilt_step
        while angle <= max_tilt + 1e-6:
            tilts.append(angle)
            angle += tilt_step
        angle = -tilt_step
        while angle >= min_tilt - 1e-6:
            tilts.append(angle)
            angle -= tilt_step
        return tilts
    if scheme == "continuous":
        tilts = []
        angle = min_tilt
        while angle <= max_tilt + 1e-6:
            tilts.append(angle)
            angle += tilt_step
        return tilts
    sem.OKBox(f"ERROR: Unknown tiltScheme '{scheme}'. Use dose_symmetric, bidirectional, or continuous.")
    sem.Exit()

def refine_geo_after_first_tilt():
    global geo, position, startTilt
    if len(geo[2]) >= 3:
        log("Refining geometry...")
        log(f"{len(geo[2])} usable CtfFind results found.")
        if len(geo[2]) >= parabolTh:
            log("Fitting paraboloid...")
            geoF = geoPara
        else:
            log("Fitting plane...")
            geoF = geoPlane
        p, cov = optimize.curve_fit(geoF, [geo[0], geo[1]], [z - geo[2][0] for z in geo[2]])
        ss = 0
        for i in range(0, len(geo[2])):
            ss += (geo[2][i] - geo[2][0] - geoF([geo[0][i], geo[1][i]], *p))**2
        rmse = np.sqrt(ss / len(geo[2]))
        log("Fit parameters: " + " # ".join(p.astype(str)))
        log(f"RMSE: {round(rmse, 3)}")
        for pos in range(0, len(position)):
            zs = geoF([position[pos][1]["SSX"], position[pos][1]["SSY"]], *p)
            z0_ref = position[pos][1]["z0"] + zs * np.cos(np.radians(startTilt)) + position[pos][1]["SSY"] * np.sin(np.radians(startTilt))
            position[pos][1]["z0"] = z0_ref
            position[pos][2]["z0"] = z0_ref
    else:
        log(f"WARNING: Not enough reliable CtfFind results ({len(geo[2])}) to refine geometry. Continuing with initial geometry model.")

def run_tilt_series(start_index):
    global tiltStepCounter
    for schedule_idx in range(start_index, total_tilt_steps):
        tilt = tilt_schedule[schedule_idx]
        tiltStepCounter = schedule_idx + 1
        log(f"\nTilt step {tiltStepCounter} out of {total_tilt_steps} ({tilt} deg)...", style=1)
        sem.SetStatusLine(1, f"Tilt step: {tiltStepCounter} / {total_tilt_steps}")
        Tilt(tilt)
        if schedule_idx == 0 and refineGeo and not recover:
            refine_geo_after_first_tilt()
        if coldFEG and tiltScheme == "dose_symmetric" and schedule_idx > 0 and schedule_idx % (2 * groupSize) == 0:
            checkColdFEG()

def Tilt(tilt):
    def calcSSChange(x, z0):                                                                    # x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
        return x[1] * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) - z0 * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

    def calcFocusChange(x, z0):                                                                 # x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
        return z0 * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) + x[1] * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

    def setTrack():
        global trackMag, origMag
        if trackDefocus < maxDefocus:
            sem.SetDefocus(position[0][pn]["focus"] + trackDefocus - targetDefocus)
        if trackExpTime > 0:
            if tilt == startTilt:
                sem.SetExposure("R", max(trackExpTime, zeroExpTime))
            else:
                sem.SetExposure("R", trackExpTime)
        if trackMag > 0:
            if tilt == startTilt:
                origMag, *_ = sem.ReportMag()
                sem.UpdateLowDoseParams("R")
            attempt = 0
            while sem.ReportMag()[0] == origMag:                                                # has to be checked, because Rec is sometimes not updated (JEOL)
                if attempt >= 10:
                    log("WARNING: Magnification could not be changed. Continuing with the same magnification for all tilt series.")
                    trackMag = 0
                    break
                sem.SetMag(trackMag)
                sem.GoToLowDoseArea("R")
                attempt += 1
            sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
            if not recover:
                sem.ImageShiftByMicrons(0, SSchange)    

    def resetTrack():
        if trackExpTime > 0:
            sem.RestoreCameraSet("R")
        if trackMag > 0:
            while sem.ReportMag()[0] != origMag:                                                # has to be checked, because Rec is sometimes not updated (JEOL)
                sem.SetMag(origMag)
                sem.GoToLowDoseArea("R")

    global recover, posResumed, resumePN, resumePlus, resumeMinus
    global is2ssMatrix, ss2isMatrix, camX, camY, c2ssMatrix, ss2cMatrix, geo
    global maxProgress, resumePercent, startTime

    # Tilt if within tilt range, skip branch if not
    if -tiltLimit <= tilt <= tiltLimit:
        skip_branch = False
        if not fixedStageTilt and tilt != startTilt:
            sem.TiltTo(tilt)
    else:
        log(f"WARNING: Tilt angle [{tilt} degrees] could not be reached. This branch of the tilt series will be aborted.")
        skip_branch = True

    if tilt < startTilt:
        increment = -step
        if not fixedStageTilt and tilt - step >= -tiltLimit:
            sem.TiltBy(-step)
            sem.TiltTo(tilt)
        pn = 2
    else:
        if not fixedStageTilt and slowTilt and startTilt < tilt <= tiltLimit:
            sem.TiltBy(-step)
            sem.TiltTo(tilt)
        increment = step
        pn = 1

    # After branch was determined, set branch to be skipped and return
    if skip_branch:
        for pos in range(len(position)):
            position[pos][pn]["skip"] = True
        return

    if fixedStageTilt:
        ensure_fixed_stage_tilt(quiet=True)

    sem.Delay(delayTilt, "s")
    realTilt = float(tilt) if fixedStageTilt else float(sem.ReportTiltAngle())

    if zeroExpTime > 0 and tilt == startTilt:
        sem.SetExposure("R", zeroExpTime)

    if recover:
        # preview align to last tracking TS
        openOldFile(targets[0]["tsfile"])
        sem.ReadFile(position[0][pn]["sec"], "O")                                               # read last image of position for AlignTo
        sem.SetDefocus(position[0][pn]["focus"])
        sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
        SSchange = 0                                                                            # needs to be defined for setTrack
        setTrack()
        if checkDewar: checkFilling()
        is_x, is_y, *_ = sem.ReportImageShift()
        sem.GoToLowDoseArea("R")
        sem.SetImageShift(0, 0)
        sem.SetImageShift(is_x, is_y)
        if beamTiltComp:
            doLafis(is_x,is_y)
        ronchi_before_preview_align("recovery preview alignment to tracking TS", pos=0, pn=pn)
        sem.L()
        if beamTiltComp:
            restoreLafis()

        alignTo("O", debug)
        bufISX, bufISY = sem.ReportISforBufferShift()
        sem.ImageShiftByUnits(position[0][pn]["ISXali"], position[0][pn]["ISYali"])             # remove accumulated buffer shifts to calculate alignment to initial startTilt image
        position[0][pn]["ISXset"], position[0][pn]["ISYset"], *_ = sem.ReportImageShift()
        for i in range(1, len(position)):
            position[i][pn]["ISXset"] += bufISX + position[0][pn]["ISXali"]                     # apply accumulated (stage dependent) buffer shifts of tracking TS to all targets
            position[i][pn]["ISYset"] += bufISY + position[0][pn]["ISYali"]
        resetTrack()
        sem.CloseFile()

        posStart = posResumed
    else:
        posStart = 0

    for pos in range(posStart, len(position)):
        log(f"\nTarget {pos + 1} / {len(position)}:", style=1)
        sem.SetStatusLine(2, "Target: " + str(pos + 1) + " / " + str(len(position)))
        if pos != 0 and position[pos][pn]["skip"]: 
            log(f"[{pos + 1}] was skipped on this branch.")
            continue
        if tilt != startTilt:
            openOldFile(targets[pos]["tsfile"])
            sem.ReadFile(position[pos][pn]["sec"], "O")                                         # read last image of position for AlignTo
        else:
            if os.path.exists(os.path.join(curDir, targets[pos]["tsfile"])):
                # Close all files incase file to be renamed is currently open
                while sem.ReportFileNumber() > 0:
                    sem.CloseFile()
                os.replace(os.path.join(curDir, targets[pos]["tsfile"]), os.path.join(curDir, targets[pos]["tsfile"]) + "~")
                log("WARNING: Tilt series file already exists. Existing file was renamed.")
            sem.OpenNewFile(targets[pos]["tsfile"])
            if not tgtPattern and "tgtfile" in targets[pos].keys():
                if refFromPreview:
                    temp_ref = os.path.splitext(targets[pos]["tgtfile"])[0] + "_tempref.mrc"
                    if os.path.exists(temp_ref) and "prevASX" in targets[pos].keys():
                        sem.ReadOtherFile(0, "O", temp_ref)                                     # reads temp reference from previewAli instead
                # BUG ? The next line overwrite what O buffer
                sem.ReadOtherFile(0, "O", targets[pos]["tgtfile"])                              # reads tgt file for first AlignTo instead

        sem.AreaForCumulRecordDose(pos + 1)                                                     # set area to accumulate record dose (counting from 1)

### Calculate and apply predicted shifts
        SSchange = 0                                                                            # only apply changes if not startTilt
        focuschange = 0
        if tilt != startTilt:
            SSchange = calcSSChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])
            focuschange = calcFocusChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

        SSYprev = position[pos][pn]["SSY"]
        SSYpred = position[pos][pn]["SSY"] + SSchange

        focuscorrection = focusSlope * (tilt - startTilt)
        position[pos][pn]["focus"] += focuscorrection
        position[pos][pn]["focus"] -= focuschange

        sem.SetDefocus(position[pos][pn]["focus"])
        if zeroDefocus != 0 and tilt == startTilt:
            sem.ChangeFocus(zeroDefocus - maxDefocus)

        sem.SetImageShift(position[pos][pn]["ISXset"], position[pos][pn]["ISYset"])
        sem.ImageShiftByMicrons(0, SSchange)

        # Apply sliding offset along tilt axis throughout tilt series when collecting montage tilt series
        if tgtMontage and (tgtTrackMnt or pos != 0) and tgtMntXOffset > 0:
            mont_offset = np.array([tgtMntXOffset * 2 * (tilt - startTilt) / (maxTilt - minTilt), 0])
            sem.ImageShiftByMicrons(mont_offset)
        else:
            mont_offset = None

### Autofocus (optional) and tracking TS settings
        if pos == 0:
            if addAF and tiltScheme == "dose_symmetric" and (tilt - startTilt) % (groupSize * step) == step and abs(tilt - startTilt) > step:
                defocus, _ = measure_defocus()
                focuserror = float(defocus) - targetDefocus
                for i in range(0, len(position)):
                    position[i][pn]["focus"] -= focuserror
                sem.SetDefocus(position[pos][pn]["focus"])

            setTrack()

### Record
        if checkDewar: checkFilling()
        checkValves()
        sem.SetFrameNameFormat(0, 0, 0x40)                                                      # turn off Sequential number
        sem.SetFrameNameFormat(0, 1, 0x400)                                                     # turn on tilt angle
        sem.SetFrameBaseName(0, 1, 0, os.path.splitext(targets[pos]["tsfile"])[0] + f"_tilt_{str(tiltStepCounter).zfill(3)}_angle")  # include collection order and tilt angle in frame name
        if beamTiltComp:
            is_x, is_y, *_ = sem.ReportImageShift()
            doLafis(is_x,is_y)
        recordWithRonchi(
            set_track_fn=setTrack if pos == 0 else None,
            acquire_label=f"Tilt {tilt:.1f} deg step {tiltStepCounter} target {pos + 1}/{len(targets)}",
            pos=pos,
            pn=pn,
        )

        bufISXpre = 0                                                                           # only non 0 if two tracking images are taken
        bufISYpre = 0
        if tilt != startTilt or (not tgtPattern and "tgtfile" in targets[pos].keys() and not noZeroRecAli): # align to previous image if it exists 
            if pos != 0: 
                sem.LimitNextAutoAlign(alignLimit)                                              # gives maximum distance for AlignTo to avoid runaway tracking
            alignTo("O", debug)
            if trackTwice and pos == 0:                                                         # track twice if alignLimit for tracking area is surpassed
                ASX, ASY = sem.ReportAlignShift()[4:6]
                if abs(ASX) > alignLimit * 1000 or abs(ASY) > alignLimit * 1000:
                    bufISXpre, bufISYpre = sem.ReportISforBufferShift()                         # have to be added only to ISset but not ISali (since ali only considers the IS chain of ali images)
                    recordWithRonchi(
                        set_track_fn=setTrack if pos == 0 else None,
                        acquire_label=f"Tilt {tilt:.1f} deg step {tiltStepCounter} target 0 (re-track)",
                        pos=pos,
                        pn=pn,
                    )
                    alignTo("O", debug)

        bufISX, bufISY = sem.ReportISforBufferShift()

        # Subtract montage offset if given
        if mont_offset is not None:
            bufISX, bufISY = np.array([bufISX, bufISY]) - ss2isMatrix @ mont_offset

        sem.ImageShiftByUnits(position[pos][pn]["ISXali"], position[pos][pn]["ISYali"])         # remove accumulated buffer shifts to calculate alignment to initial startTilt image

        if beamTiltComp:
            restoreLafis()

        position[pos][pn]["ISXset"], position[pos][pn]["ISYset"], *_ = sem.ReportImageShift()
        position[pos][pn]["SSX"], position[pos][pn]["SSY"] = sem.ReportSpecimenShift()

        # Collect surrounding tiles for montage tilt series
        if tgtMontage and (tgtTrackMnt or pos != 0):
            sem.ImageShiftByUnits(-bufISX - position[pos][pn]["ISXali"], -bufISY - position[pos][pn]["ISYali"]) # reset shifts to already taken center image
            for i in range(-tgtMntSize, tgtMntSize + 1):
                for j in range(-tgtMntSize, tgtMntSize + 1):
                    if i == j == 0: continue
                    if tilt != startTilt:
                        openOldFile(os.path.splitext(targets[pos]["tsfile"])[0] + "_" + str(i) + "_" + str(j) + ".mrc")
                    else:
                        sem.OpenNewFile(os.path.splitext(targets[pos]["tsfile"])[0] + "_" + str(i) + "_" + str(j) + ".mrc")

                    montX, montY = (i - i * tgtMntOverlap) * min([camX, camY]), (j - j * tgtMntOverlap) * min([camX, camY])
                    pixelShiftFromCenter = f"{montX} {montY}"

                    # Apply sliding offset along tilt axis throughout tilt series
                    if tgtMntXOffset > 0:
                        # Adjust x shift by tgtMntXOffset [microns at max tilt] by fraction of tilt series along SSX
                        montX, montY = np.array([montX, montY]) + ss2cMatrix @ np.array([tgtMntXOffset * 2 * (tilt - startTilt) / (maxTilt - minTilt), 0])

                    sem.ImageShiftByPixels(montX, montY)
                    if tgtMntFocusCor:
                        montSSX, montSSY = c2ssMatrix @ np.array([montX, montY])

                        # With sample geometry (needs to be tested)
                        correctedFocus = position[pos][pn]["focus"] - np.cos(np.radians(realTilt)) * np.tan(np.radians(pretilt)) * (np.cos(np.radians(rotation)) / np.cos(np.radians(realTilt)) * montSSY - np.sin(np.radians(rotation)) * montSSX) - np.tan(np.radians(realTilt)) * montSSY 
                        # Without sample geometry
                        #correctedFocus = position[pos][pn]["focus"] - np.tan(np.radians(realTilt)) * montSSY

                    sem.SetDefocus(correctedFocus)
                    if beamTiltComp:
                        is_x,is_y,*_ = sem.ReportImageShift()
                        doLafis(is_x,is_y)
                    recordWithRonchi(
                        set_track_fn=setTrack if pos == 0 else None,
                        run_ronchi=ronchiMontage,
                        acquire_label=f"Montage tile ({i},{j}) tilt {tilt:.1f} deg target {pos + 1}/{len(targets)}",
                        pos=pos,
                        pn=pn,
                    )

                    mont_SSX, mont_SSY = sem.ReportSpecimenShift()

                    sem.ImageShiftByPixels(-montX, -montY)
                    if beamTiltComp:
                        restoreLafis()

                    # Add shift to all montage tilt series mdoc files for auto stitching
                    if extendedMdoc:
                        sem.AddToAutodoc("PixelShiftFromCenter", pixelShiftFromCenter)
                        sem.WriteAutodoc()

                    sem.CloseFile()

        position[pos][pn]["focus"] -= focuscorrection                                           # remove correction or it accumulates

        dose = sem.ImageConditions("A")[0]
        if dose > 0:
            sem.AccumulateRecordDose(dose)
            position[pos][1]["dose"] += dose
            position[pos][2]["dose"] += dose

        if pos == 0:                                                                            # apply measured shifts of first/tracking position to other positions
            for i in range(1, len(position)):
                position[i][pn]["ISXset"] += bufISX + bufISXpre + position[pos][pn]["ISXali"]   # apply accumulated (stage dependent) buffer shifts of tracking TS to all targets
                position[i][pn]["ISYset"] += bufISY + bufISYpre + position[pos][pn]["ISYali"]
                if tilt == startTilt:                                                           # also save shifts from startTilt image for second branch since it will alignTo the startTilt image
                    position[i][2]["ISXset"] += bufISX + bufISXpre
                    position[i][2]["ISYset"] += bufISY + bufISYpre
            if tilt == startTilt:                                                               # do not forget about 0 position
                position[0][2]["ISXset"] += bufISX + bufISXpre
                position[0][2]["ISYset"] += bufISY + bufISYpre

            resetTrack()

        position[pos][pn]["ISXali"] += bufISX
        position[pos][pn]["ISYali"] += bufISY
        if tilt == startTilt:                                                                   # save alignment of first tilt to tgt file for the second branch
            position[pos][2]["ISXali"] += bufISX
            position[pos][2]["ISYali"] += bufISY

        aErrX, aErrY = is2ssMatrix @ np.array([position[pos][pn]["ISXali"], position[pos][pn]["ISYali"]])

        log(f"[{pos + 1}] Prediction: y = {round(SSYpred, 3)} microns | z = {round(position[pos][pn]['focus'], 3)} microns | z0 = {round(position[pos][pn]['z0'], 3)} microns")
        log(f"[{pos + 1}] Reality: y = {round(position[pos][pn]['SSY'], 3)} microns")
        log(f"[{pos + 1}] Focus change: {round(focuschange, 3)} microns | Focus correction: {round(focuscorrection, 3)} microns")
        log(f"[{pos + 1}] Alignment error: x = {round(aErrX * 1000)} nm | y = {round(aErrY * 1000)} nm")        

### Calculate new z0

        ddy = position[pos][pn]["SSY"] - SSYprev
        if (tilt == startTilt or
                (ignoreNegStart and pn == 2 and len(position[pos][pn]["shifts"]) == 0) or
                recover or
                (resumePN == 1 and tilt == resumePlus + step and pos < posResumed) or
                (resumePN == 1 and tilt == resumeMinus - step) or
                (resumePN == 2 and tilt == resumeMinus - step and pos < posResumed) or
                (resumePN == 2 and tilt == resumePlus + step)):        
                # ignore shift if first image or first shift of second branch or first image after resuming run (all possible conditions)
            ddy = calcSSChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

        position[pos][pn]["shifts"].append(ddy)
        position[pos][pn]["angles"].append(realTilt)

        if len(position[pos][pn]["shifts"]) > dataPoints:
            position[pos][pn]["shifts"].pop(0)
            position[pos][pn]["angles"].pop(0)

        position[pos][pn]["z0"], cov = optimize.curve_fit(calcSSChange, np.vstack((position[pos][pn]["angles"], [position[pos][pn]["n0"] for i in range(0, len(position[pos][pn]["angles"]))])), position[pos][pn]["shifts"], p0=(position[pos][pn]["z0"]))
        position[pos][pn]["z0"] = position[pos][pn]["z0"][0]

        if doCtfFind:
            try:
                sem.NoMessageBoxOnError(1)
                cfind = sem.CtfFind("A", (min(maxDefocus, trackDefocus) - 2), min(-0.2, minDefocus + 2))
                log(f"[{pos + 1}] CtfFind: {round(cfind[0], 3)} microns ({round(cfind[-1], 2)} A)")
            except:
                log(f"WARNING: CtfFind crashed on {targets[pos]['tsfile']} section {int(sem.ReportFileZsize()) - 1}. Trying to continue...")
            finally:
                sem.NoMessageBoxOnError(0)

        if doCtfPlotter:
            try:
                sem.NoMessageBoxOnError(1)
                cplot = sem.Ctfplotter("A", (min(maxDefocus, trackDefocus) - 2), min(-0.2, minDefocus + 2), 1, 0, pretilt)
                log(f"[{pos + 1}] Ctfplotter: {round(cplot[0], 3)} microns")
            except:
                log(f"WARNING: Ctfplotter crashed on {targets[pos]['tsfile']} section {int(sem.ReportFileZsize()) - 1}. Trying to continue...")
            finally:
                sem.NoMessageBoxOnError(0)

        if refineGeo and tilt == startTilt:
            if doCtfPlotter:
                geo[0].append(position[pos][pn]["SSX"])
                geo[1].append(position[pos][pn]["SSY"])
                geo[2].append(cplot[0])
            elif doCtfFind and len(cfind) > 5:
                if cfind[5] < fitLimit:                                                         # save vectors for refineGeo only if CTF fit has reasonable resolution
                    geo[0].append(position[pos][pn]["SSX"])
                    geo[1].append(position[pos][pn]["SSY"])
                    geo[2].append(cfind[0])

        position[pos][pn]["sec"] = int(sem.ReportFileZsize()) - 1                               # save section number for next alignment

        # progress = collected images * (positions - skipped positions) + current position - skipped positions scaled assuming homogeneous distribution of skipped positions
        progress = position[pos][pn]["sec"] * (len(position) - skippedTgts) + pos - skippedTgts * pos / len(position) + 1
        percent = round(100 * (progress / maxProgress), 1)
        bar = '#' * int(percent / 2) + '_' * (50 - int(percent / 2))
        if percent - resumePercent > 0:
            remTime = int((sem.ReportClock() - startTime - dewarFillTime) / (percent - resumePercent) * (100 - percent) / 60)
        else:
            remTime = "?"
        log(f"Progress: |{bar}| {percent} % ({remTime} min remaining)")

        if extendedMdoc:
            sem.AddToAutodoc("SpecimenShift", str(position[pos][pn]["SSX"]) + " " + str(position[pos][pn]["SSY"]))
            sem.AddToAutodoc("EucentricOffset", str(position[pos][pn]["z0"]))
            if tgtMontage:
                sem.AddToAutodoc("PixelShiftFromCenter", "0 0")
            if doCtfFind:
                sem.AddToAutodoc("CtfFind", str(cfind[0]))
            if doCtfPlotter:
                sem.AddToAutodoc("Ctfplotter", str(cplot[0]))
            sem.WriteAutodoc()

        sem.CloseFile()

### Abort conditions
        if np.linalg.norm(np.array([position[pos][pn]["SSX"], position[pos][pn]["SSY"]], dtype=float)) > imageShiftLimit - alignLimit:
            position[pos][pn]["skip"] = True
            log(f"WARNING: Target [{pos + 1}] is approaching the image shift limit. This branch will be aborted.")

        if minCounts > 0:
            meanCounts = sem.ReportMeanCounts()
            expTime, *_ = sem.ReportExposure("R")
            if meanCounts / expTime < minCounts:
                position[pos][pn]["skip"] = True
                log(f"WARNING: Target [{pos + 1}] was too dark. This branch will be aborted.")

        if tilt >= maxTilt or tilt <= minTilt:
            position[pos][pn]["skip"] = True
            if maxTilt - startTilt != abs(minTilt - startTilt):
                log(f"WARNING: Target [{pos + 1}] has reached the final tilt angle. This branch will be aborted.")            

        updateTargets(runFileName, targets, position, position[pos][pn]["sec"], pos)    

### Refine energy filter slit if appropriate
    if tgtPattern and slitInterval > 0 and (lastSlitCheck - sem.ReportClock() / 60) > slitInterval:
        checkSlit(np.array([vecB0, vecB1]), size, realTilt, pn)

    if zeroExpTime > 0 and tilt == startTilt:
        sem.RestoreCameraSet("R")

    if recover:
        recover = False    



def find_tgts_files(fileStem, curDir):
    tf = sorted(glob.glob(os.path.join(curDir, fileStem + ".txt")))
    tfr = sorted(glob.glob(os.path.join(curDir, fileStem + "_run??.txt")))
    tf.extend(tfr)
    return tf


def resolve_nav_items():
    if nav_item_list:
        return [int(x) for x in nav_item_list]
    sem.ReportNavItem()
    return [int(sem.GetVariable("navIndex"))]


def validate_nav_lists():
    if not nav_item_list:
        return
    n = len(nav_item_list)
    if nav_pretilt_list and len(nav_pretilt_list) != n:
        sem.OKBox(f"ERROR: nav_pretilt_list length ({len(nav_pretilt_list)}) must match nav_item_list ({n}).")
        sem.Exit()
    if nav_rotation_list and len(nav_rotation_list) != n:
        sem.OKBox(f"ERROR: nav_rotation_list length ({len(nav_rotation_list)}) must match nav_item_list ({n}).")
        sem.Exit()
    if nav_start_defocus_list and len(nav_start_defocus_list) != n:
        sem.OKBox(f"ERROR: nav_start_defocus_list length ({len(nav_start_defocus_list)}) must match nav_item_list ({n}).")
        sem.Exit()


def get_nav_start_defocus(item_index):
    """Objective defocus [um] to set at nav item start (before autofocus)."""
    if nav_start_defocus_list and item_index < len(nav_start_defocus_list):
        return float(nav_start_defocus_list[item_index])
    return float(default_nav_start_defocus)


def apply_nav_start_defocus(item_index, when=""):
    """Set objective defocus before autofocus so measure/autofocus is not started too far off."""
    defocus = get_nav_start_defocus(item_index)
    sem.SetDefocus(defocus)
    suffix = f" {when}" if when else ""
    log(f"NOTE: Set objective defocus to {defocus:.2f} um at nav item start{suffix}.")


def get_nav_geometry(item_index):
    if nav_item_list:
        pt = nav_pretilt_list[item_index] if item_index < len(nav_pretilt_list) else default_pretilt
        rot = nav_rotation_list[item_index] if item_index < len(nav_rotation_list) else default_rotation
        return float(pt), float(rot)
    return default_pretilt, default_rotation


def nav_run_status(tf_path):
    if tf_path is None:
        return "fresh"
    basename = os.path.basename(tf_path)
    if "_run" not in basename:
        return "fresh"
    resume_sec, resume_pos = 0, 0
    with open(tf_path) as f:
        for ln in f:
            if ln.startswith("_spos"):
                col = ln.strip().split(" ")
                parts = col[2].split(",")
                resume_sec = int(parts[0])
                resume_pos = int(parts[1])
                break
    if resume_sec > 0 or resume_pos > 0:
        return "incomplete"
    return "complete"


def get_nav_tgts_path(nav_idx, curDir=None):
    sem.SetSelectedNavItem(nav_idx)
    sem.ReportNavItem()
    navNote = sem.GetVariable("navNote")
    fileStem, fileExt = os.path.splitext(navNote)
    if fileStem == "" or fileExt != ".txt":
        return None
    if curDir is None:
        curDir = sem.ReportDirectory()
    tf = find_tgts_files(fileStem, curDir)
    if not tf:
        return None
    return tf[-1]


def batch_resume_start_index(items_to_run):
    curDir = sem.ReportDirectory()
    for i, nav_idx in enumerate(items_to_run):
        tf_path = get_nav_tgts_path(nav_idx, curDir)
        if tf_path is None:
            return i
        status = nav_run_status(tf_path)
        if status in ("incomplete", "fresh"):
            return i
    return len(items_to_run)


def _reset_ronchi_to_session_start(when=""):
    """Restore XLens(2) and C3 to values captured at batch start in checkRonchigramSetup."""
    if not doRonchigram:
        return
    suffix = f" {when}" if when else ""
    if hasXLens and ronchiStartXLensX is not None and ronchiStartXLensY is not None:
        sem.SetXLensDeflector(2, ronchiStartXLensX, ronchiStartXLensY)
        log(
            f"NOTE: Reset XLensDeflector(2) to session start "
            f"({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}){suffix}."
        )
    if ronchiStartC3Offset is not None:
        sem.SetImageDistanceOffset(ronchiStartC3Offset)
        log(f"NOTE: Reset ImageDistanceOffset to session start {ronchiStartC3Offset:.2f} um{suffix}.")


def reset_scope_between_nav_items():
    """Reset scope state between multi-nav items so a bad area does not carry over."""
    while sem.ReportFileNumber() > 0:
        sem.CloseFile()
    sem.TiltTo(0)
    sem.SetImageShift(0, 0)
    if trackMag > 0:
        try:
            sem.RestoreLowDoseParams("R")
        except (AttributeError, TypeError):
            pass
    _reset_ronchi_to_session_start("between nav items")


def run_one_nav_item(nav_idx, item_index, batch_recover=False, batch_recover_accepted=False):
    """Run full PACEtomo acquisition for one Navigator item."""
    global startTilt, minTilt, maxTilt, pretilt, rotation, recover, recoverInput, realign
    global targets, savedRun, resume, geoPoints, position, navID, navNote, fileStem, curDir
    global runFileName, targetDefocus, tf, skippedTgts, focus0, tiltLimit, vecA0, vecA1, vecB0, vecB1, size
    global minDefocus, maxDefocus, branchsteps, tilt_schedule, total_tilt_steps, schedule_start
    global geo, posResumed, resumePN, resumePlus, resumeMinus, resumePercent, maxProgress, startTime
    global tiltStepCounter, dewarFillTime, lastSlitCheck
    global is2ssMatrix, ss2isMatrix, s2ssMatrix, camX, camY, c2ssMatrix, ss2cMatrix, origMag
    global positionFocus, minFocus0, ISX0, ISY0, SSX0, SSY0

    startTilt = default_startTilt
    minTilt = default_minTilt
    maxTilt = default_maxTilt
    pretilt = default_pretilt
    rotation = default_rotation
    recover = False
    recoverInput = 0
    realign = False

    ### Find target file
    sem.SetSelectedNavItem(nav_idx)
    sem.ReportNavItem()
    navID = int(sem.GetVariable("navIndex"))
    navNote = sem.GetVariable("navNote")
    fileStem, fileExt = os.path.splitext(navNote)
    curDir = sem.ReportDirectory()

    if fileStem != "" and fileExt == ".txt":
        tf = find_tgts_files(fileStem, curDir)
        while tf == []:
            searchInput = sem.YesNoBox("\n".join(["Target file not found! Please choose the directory containing the target file!", "WARNING: All future target files will be searched here!"]))
            if searchInput == 0:
                return False
            sem.UserSetDirectory("Please choose the directory containing the target file!")
            curDir = sem.ReportDirectory()
            tf = find_tgts_files(fileStem, curDir)
    else:
        log(f"ERROR: Nav item {nav_idx} note does not contain a target file.")
        return False

    # Check if frame folder is set reasonably
    if sem.ReportCameraProperty(0, "K2Type") > 0:
        framePath = sem.ReportFrameSavingPath()
        framePar = os.path.abspath(os.path.join(framePath, os.pardir))
        curPar = os.path.abspath(os.path.join(curDir, os.pardir))
        if (framePath == "NONE" or (framePath not in curDir and curDir not in framePath and framePar not in curDir and curPar not in framePath)) and sem.IsVariableDefined("warningFramePath") == 0:
            sem.Pause("WARNING: Current frame path (" + framePath + ") does not seem plausible or camera does not save frames.")
            sem.SetPersistentVar("warningFramePath", "")
    else:
        log("WARNING: Camera frame path could not be obtained for your camera.")

    sem.SaveLogOpenNew(navNote.split("_tgts")[0])

    log(f"PACEtomo Version {versionPACE}", color=5, style=1)
    sem.ProgramTimeStamps()

    # Open last tgts or tgts_run file and read contents
    targets, savedRun, resume, geoPoints = parseTargets(tf[-1])
    pretilt, rotation = get_nav_geometry(item_index)


    # Sanity check of settings
    if maxDefocus > minDefocus:
        minDefocus, maxDefocus = maxDefocus, minDefocus
    if maxTilt < minTilt:
        minTilt, maxTilt = maxTilt, minTilt
    if (maxTilt - startTilt) % step != 0 or (startTilt - minTilt) % step != 0:
        maxTilt = round(int((maxTilt - startTilt) / step) * step + startTilt, 1)
        minTilt = round(int((minTilt - startTilt) / step) * step + startTilt, 1)
        print(f"WARNING: Tilt increment does not divide evenly into tilt range. Tilt range will be adjusted to: {minTilt}, {maxTilt}")

    ### Recovery data
    recoverInput = 0
    recover = False
    realign = False
    if savedRun != False and (resume["sec"] > 0 or resume["pos"] > 0):
        if batch_recover:
            recoverInput = 1 if batch_recover_accepted else 0
        elif nav_item_list:
            recoverInput = 0
            sem.AllowFileOverwrite(1)
        else:
            recoverInput = sem.YesNoBox("The target file contains recovery data. Do you want to attempt to continue the acquisition? Tracking accuracy might be impacted.")
        if recoverInput == 1:
            recover = True
            while sem.ReportFileNumber() > 0:
                sem.CloseFile()

            stageX, stageY, stageZ = sem.ReportStageXYZ()
            if abs(stageX - float(targets[0]["stageX"])) > 1.0 or abs(stageY - float(targets[0]["stageY"])) > 1.0: # test if stage was moved (with 1 micron wiggle room)
                userRealign = sem.YesNoBox("It seems that the stage was moved since stopping acquisition. Do you want to realign to the tracking target before resuming? This will also reset prediction parameters reducing tracking accuracy.")    
                realign = True if userRealign == 1 else False
        else:
            sem.AllowFileOverwrite(1)

    ### Start setup
    dumpVars(os.path.splitext(os.path.basename(tf[-1]))[0])                                         # write settings vars to text file

    sem.ResetClock()

    targetDefocus = maxDefocus                                                                      # use highest defocus for tracking TS
    sem.SetTargetDefocus(targetDefocus)

    # Collect exposure settings
    expTime = sem.ReportExposure("R")[0]

    if recover:
        log("##### Recovery attempt of PACEtomo with parameters: #####", style=1)
    else:
        log("##### Starting new PACEtomo with parameters: #####", style=1)
    log(f"Start: {startTilt} deg - Min/Max: {minTilt}/{maxTilt} deg ({step} deg increments)")
    log(f"Tilt scheme: {tiltScheme}" + (f" (groupSize={groupSize})" if tiltScheme == "dose_symmetric" else ""))
    log(f"Data points used: {dataPoints}")
    log(f"Target defocus range (min/max/step): {minDefocus}/{maxDefocus}/{stepDefocus}")
    log(f"Defocus method (measure / autofocus): {defocusMethod}")
    if fixedStageTilt:
        log(
            f"fixedStageTilt: stage held at {fixedStageTiltAngle:.1f} deg; "
            f"scheduled series {minTilt} to {maxTilt} deg ({step} deg step)"
        )
    log(f"Sample pretilt (rotation): {pretilt} ({rotation})")
    log(f"Nav item start defocus: {get_nav_start_defocus(item_index):.2f} um")
    log(f"Tilt axis offset: {round(tiltAxisOffset, 3)}")
    log(f"Focus correction slope: {focusSlope}")
    log(f"Exposure time per tilt: {round(expTime, 3)} s (total: {round(expTime * int((maxTilt - minTilt) / step), 3)} s)")

    if trackMag > 0:
        log("WARNING: A magnification offset for the tracking target changes the Low Dose Record mode temporarily. Please double-check your Record mode in case the script is stopped prematurely or crashes!")

    if startTilt * pretilt > 0:
        log("WARNING: Start tilt and pretilt have the same sign! If you want to compensate for the pretilt, the start tilt should have the opposite sign!")

    ### Create run file
    counter = 1
    while os.path.exists(os.path.join(curDir, fileStem + "_run" + str(counter).zfill(2) + ".txt")):
        counter += 1
    runFileName = os.path.join(curDir, fileStem + "_run" + str(counter).zfill(2) + ".txt")

    ### Initial actions
    if not recover:
        log("Moving to target area...")

        sem.SetCameraArea("V", "F")                                                                 # set View to Full for Eucentricity
        sem.MoveToNavItem(navID)
        apply_nav_start_defocus(item_index, "after MoveToNavItem")
        log("Refining eucentricity...")
        sem.GoToLowDoseArea("V")
        is_x, is_y, *_ = sem.ReportImageShift()
        sem.SetImageShift(0, 0)
        sem.Eucentricity(1)
        sem.SetImageShift(is_x, is_y)
        sem.UpdateItemZ()
        sem.RestoreCameraSet("V")

        log("Realigning to target 1...")
        if alignToP:
            x, y, binning, exp, *_ = sem.ImageProperties("P")
            sem.SetExposure("V", exp)
            sem.SetBinning("V", int(binning))
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            sem.CropCenterToSize("A", int(x), int(y))
            alignTo("P", debug)
            sem.RestoreCameraSet("V")
            if refineVec and tgtPattern and size is not None:
                if float(sem.ReportDefocus()) < -50:
                    log("WARNING: Large defocus offsets for View can cause additional offsets in image shift upon mag change.")
                size = int(size)
                log("Refining target pattern...")
                is_x, is_y, *_ = sem.ReportImageShift()
                sem.GoToLowDoseArea("R")
                sem.SetImageShift(0, 0)
                sem.SetImageShift(is_x, is_y)
                ISX0, ISY0, *_ = sem.ReportImageShift()
                SSX0, SSY0 = sem.ReportSpecimenShift()
                log(f"Vector A: ({vecA0}, {vecA1})")
                shiftx = size * vecA0
                shifty = size * vecA1
                is_x, is_y, *_ = sem.ReportImageShift()
                sem.GoToLowDoseArea("V")
                sem.SetImageShift(0, 0)
                sem.SetImageShift(is_x, is_y)
                sem.ImageShiftByMicrons(shiftx, shifty)
                sem.V()
                alignTo("P", debug)
                is_x, is_y, *_ = sem.ReportImageShift()
                sem.GoToLowDoseArea("R")
                sem.SetImageShift(0, 0)
                sem.SetImageShift(is_x, is_y)

                SSX, SSY = sem.ReportSpecimenShift()
                SSX -= SSX0
                SSY -= SSY0        
                if np.linalg.norm([shiftx - SSX, shifty - SSY]) > 0.5:
                    log("WARNING: Refined vector differs by more than 0.5 microns! Original vectors will be used.")
                else:
                    vecA0, vecA1 = (round(SSX / size, 4), round(SSY / size, 4))
                    log(f"Refined vector A: ({vecA0}, {vecA1})")

                    sem.SetImageShift(ISX0, ISY0)                                                   # reset IS to center position
                    log(f"Vector B: ({vecB0}, {vecB1})")
                    shiftx = size * vecB0
                    shifty = size * vecB1
                    sem.ImageShiftByMicrons(shiftx, shifty)

                    is_x, is_y, *_ = sem.ReportImageShift()
                    sem.GoToLowDoseArea("V")
                    sem.SetImageShift(0, 0)
                    sem.SetImageShift(is_x, is_y)
                    sem.V()
                    alignTo("P", debug)
                    is_x, is_y, *_ = sem.ReportImageShift()
                    sem.GoToLowDoseArea("R")
                    sem.SetImageShift(0, 0)
                    sem.SetImageShift(is_x, is_y)
                    SSX, SSY = sem.ReportSpecimenShift()
                    SSX -= SSX0
                    SSY -= SSY0
                    if np.linalg.norm([shiftx - SSX, shifty - SSY]) > 0.5:
                        log("WARNING: Refined vector differs by more than 0.5 microns! Original vectors will be used.")
                    else:
                        vecB0, vecB1 = (round(SSX / size, 4), round(SSY / size, 4))
                        log(f"Refined vector B: ({vecB0}, {vecB1})")

                        targetNo = 0
                        for i in range(-size,size+1):
                            for j in range(-size,size+1):
                                if i == j == 0: continue
                                targetNo += 1
                                SSX = i * vecA0 + j * vecB0
                                SSY = i * vecA1 + j * vecB1
                                targets[targetNo]["SSX"] = str(SSX)
                                targets[targetNo]["SSY"] = str(SSY)
                        log("NOTE: Target pattern was overwritten using refined vectors.")
                sem.SetImageShift(ISX0, ISY0)                                                       # reset IS to center position
        else:
            #sem.RealignToOtherItem(navID, 1) # <= sometimes unreliable
            realignTo(nav_id=navID, target=targets[0])

        if measureGeo:
            log("Measuring geometry...")
            if fixedStageTilt:
                ensure_fixed_stage_tilt("for measureGeo")
            elif int(round(float(sem.ReportTiltAngle()))) != 0:
                sem.TiltTo(0)
            if len(geoPoints) > 0 and "SSX" in geoPoints[0].keys():                                 # if there are geo points in tgts file, adjust format from dict to list
                geoPoints = [[point["SSX"], point["SSY"]] for point in geoPoints]
            if len(geoPoints) < 3 and tgtPattern and size is not None:
                if size > 1:
                    geoPoints.append([0.5 * (vecA0 + vecB0), 0.5 * (vecA1 + vecB1)])
                geoPoints.append([(size - 0.5) * (vecA0 + vecB0), (size - 0.5) * (vecA1 + vecB1)])
                geoPoints.append([(size - 0.5) * (vecA0 - vecB0), (size - 0.5) * (vecA1 - vecB1)])
                geoPoints.append([(size - 0.5) * (-vecA0 + vecB0), (size - 0.5) * (-vecA1 + vecB1)])
                geoPoints.append([(size - 0.5) * (-vecA0 - vecB0), (size - 0.5) * (-vecA1 - vecB1)])

            # Clean geo_points beyond image shift limit
            geoPoints = [point for point in geoPoints if np.linalg.norm(np.array([point[0], point[1]], dtype=float)) < imageShiftLimit]

            if len(geoPoints) >= 3:
                geoXYZ = [[], [], []]
                sem.GoToLowDoseArea("R")
                ISX0, ISY0, *_ = sem.ReportImageShift()
                for i in range(len(geoPoints)):
                    sem.ImageShiftByMicrons(geoPoints[i][0], geoPoints[i][1])
                    defocus, drift = measure_defocus()
                    drift_ok = (
                        defocusMethod in ("ctf", "beam_tilt_sem")
                        or np.linalg.norm(drift) >= 0.01
                    )
                    if abs(defocus) >= 0.01 and drift_ok:
                        geoXYZ[0].append(geoPoints[i][0])
                        geoXYZ[1].append(geoPoints[i][1])
                        geoXYZ[2].append(defocus)
                    else:
                        log("WARNING: Measured defocus is 0. This geo point will not be considered.")
                    sem.SetImageShift(ISX0, ISY0)                                                   # reset IS to center position
                if len(geoXYZ[0]) >= 3:
                    ##########
                    # Source: https://math.stackexchange.com/q/99317
                    # subtract out the centroid and take the SVD, extract the left singular vectors, the corresponding left singular vector is the normal vector of the best-fitting plane
                    svd = np.linalg.svd(geoXYZ - np.mean(geoXYZ, axis=1, keepdims=True))
                    left = svd[0]
                    norm = left[:, -1]
                    ##########        
                    log(f"Fitted plane into cloud of {len(geoXYZ[0])} points ({len(geoPoints) - len(geoXYZ[0])} discarded).")
                    log(f"Normal vector: {norm}")

                    # Errors
                    errors = []
                    for point in zip(*geoXYZ):
                        errors.append(np.dot(norm, point - np.mean(geoXYZ, axis=1)) ** 2)
                    log(f"Fitting error: {np.mean(errors)}")

                    if debug:
                        log("DEBUG:\nGeo points [x, y, z, err]:")
                        for point in zip(*geoXYZ, errors):
                            log(f"# {point}", color=1)

                    # Calculate pretilt and rotation
                    sign = 1 if norm[1] <= 0 else -1
                    pretilt = round(sign * np.degrees(np.arccos(norm[2])), 1)
                    log(f"Estimated pretilt: {pretilt} degrees", style=1)
                    rotation = round(-np.degrees(np.arctan(norm[0]/norm[1])), 1)
                    log(f"Estimated rotation: {rotation} degrees", style=1)

                    if startTilt * pretilt > 0:
                        log("WARNING: Start tilt and pretilt have the same sign! If you want to compensate for the pretilt, the start tilt should have the opposite sign!")
                else:
                    log("WARNING: Not enough geo points could be checked successfully. Geometry could not be measured.")
            else:
                log("WARNING: Not enough geo points were defined. Geometry could not be measured.")

        if autoStartTilt or tiltTargets != 0:
            startTiltOri = startTilt
            if tiltTargets != 0:
                # Use tilt at which targets were selected as startTilt
                startTilt = tiltTargets
            elif autoStartTilt:
                # Adjust start tilt to compensate for measured pretilt
                startTilt = -int(round(np.degrees(np.arctan(np.sin(np.radians(pretilt)) * np.cos(np.radians(rotation)) / np.cos(np.radians(pretilt))))))
            maxTilt = np.clip(maxTilt - startTiltOri + startTilt, -int(tiltLimit), int(tiltLimit))
            minTilt = np.clip(minTilt - startTiltOri + startTilt, -int(tiltLimit), int(tiltLimit))

            log("WARNING: Automatically adjusted tilt series parameters!")
            log(f"Start: {startTilt} deg - Min/Max: {minTilt}/{maxTilt} deg ({step} deg increments)", style=1)

            # Update branch steps
            branchsteps = max(maxTilt - startTilt, abs(minTilt - startTilt)) / groupSize / step

        tilt_schedule = build_tilt_schedule(tiltScheme, startTilt, minTilt, maxTilt, step, groupSize)
        total_tilt_steps = len(tilt_schedule)
        if tiltScheme == "continuous" and abs(startTilt - minTilt) > 1e-6:
            log(f"NOTE: Continuous scheme collects {minTilt} to {maxTilt} deg; startTilt ({startTilt}) is not the first angle.")
        log(f"Tilt series: {total_tilt_steps} angles ({tilt_schedule[0]} to {tilt_schedule[-1]} deg)")

        first_tilt = tilt_schedule[0]
        if fixedStageTilt:
            log(
                f"fixedStageTilt: skipping tilt-to-start; "
                f"stage stays at {fixedStageTiltAngle:.1f} deg "
                f"(first scheduled angle {first_tilt:.1f} deg)"
            )
            ensure_fixed_stage_tilt("before tilt series")
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            sem.Copy("A", "O")
        else:
            log("Tilting to start tilt angle...")
            # backlash correction
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            sem.Copy("A", "O")

            curTilt = int(round(float(sem.ReportTiltAngle())))

            # Walk up if necessary
            while abs(first_tilt - curTilt) > 10:
                log(f"DEBUG: Doing walkup to {curTilt + (10 if first_tilt > curTilt else -10)}...")
                sem.TiltTo(curTilt + (10 if first_tilt > curTilt else -10))
                is_x, is_y, *_ = sem.ReportImageShift()
                sem.GoToLowDoseArea("V")
                sem.SetImageShift(0, 0)
                sem.SetImageShift(is_x, is_y)
                sem.V()
                alignTo("O", debug)

                sem.V()
                sem.Copy("A", "O")
                curTilt = int(round(float(sem.ReportTiltAngle())))

            sem.TiltTo(first_tilt - step)
            sem.TiltTo(first_tilt)

        sem.V()
        alignTo("O", debug)
        is_x, is_y, *_ = sem.ReportImageShift()
        sem.GoToLowDoseArea("R")
        sem.SetImageShift(0, 0)
        sem.SetImageShift(is_x, is_y)

        if not tgtPattern and previewAli:
            sem.SetImageShift(is_x, is_y)
            if beamTiltComp:
                doLafis(is_x,is_y)
            ronchi_before_preview_align("tracking target map preview alignment")
            sem.LoadOtherMap(navID, "O")                                                            # preview ali before first tilt image is taken
            #TODO AcquiteToMatchBuffer forces the scope params, including xt
            # to be the same as the one in the buffer (???).  If so, would restore pre-lafis
            sem.AcquireToMatchBuffer("O")                                                           # in case view image was saved for tracking target
            if beamTiltComp:
                restoreLafis()
            alignTo("O", debug)

        ISX0, ISY0, *_ = sem.ReportImageShift()
        SSX0, SSY0 = sem.ReportSpecimenShift()

        autofocus_apply(targetDefocus)
        focus0 = float(sem.ReportDefocus())
        positionFocus = focus0                                                                      # set maxDefocus as focus0 and add focus steps in loop
        minFocus0 = focus0 - maxDefocus + minDefocus

        is_x, is_y, *_ = sem.ReportImageShift()
        sem.GoToLowDoseArea("R")
        sem.SetImageShift(0, 0)
        sem.SetImageShift(is_x, is_y)
        s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
        is2ssMatrix = np.array(sem.ISToSpecimenMatrix(0)).reshape((2, 2))
        ss2isMatrix = np.array(sem.SpecimenToISMatrix(0)).reshape((2, 2))
        camX, camY, *_ = sem.CameraProperties()
        c2ssMatrix = np.array(sem.CameraToSpecimenMatrix(0)).reshape((2, 2))
        ss2cMatrix = np.array(sem.SpecimenToCameraMatrix(0)).reshape((2, 2))
        if debug:
            log("DEBUG: Conversion matrices:")
            log(f"    Stage to Specimen: {s2ssMatrix}", color=1)
            log(f"    IS to Specimen: {is2ssMatrix}", color=1)
            log(f"    Image to Specimen: {c2ssMatrix}", color=1)
            log(f"    Specimen to Camera: {ss2cMatrix}", color=1)

        if previewAli:
            sem.SetDefocus(min(focus0, focus0 - 5 - targetDefocus))                                 # set defocus for Preview to at least -5 micron
    ### Target setup
        log(f"Setting up {len(targets)} targets...")

        posTemplate = {"SSX": 0, "SSY": 0, "focus": 0, "z0": 0, "n0": 0, "shifts": [], "angles": [], "ISXset": 0, "ISYset": 0, "ISXali": 0, "ISYali": 0, "dose": 0, "sec": 0, "skip": False, "c3_offset": None}
        position = []
        skippedTgts = 0
        for i, tgt in enumerate(targets):
            position.append([])
            position[-1].append(copy.deepcopy(posTemplate))

            log(f"Target {i + 1}...")
            skip = False
            if "skip" in tgt.keys() and tgt["skip"] == "True":
                log(f"WARNING: Target [{str(i + 1).zfill(3)}] was set to be skipped.")
                skip = True
            if "SSX" not in tgt.keys() and "stageX" in tgt.keys():                                  # if SS coords are missing but stage coords are present, calc SS coords
                tgt["SSX"], tgt["SSY"] = s2ssMatrix @ np.array([float(tgt["stageX"]) - float(targets[0]["stageX"]), float(tgt["stageY"]) - float(targets[0]["stageY"])])
            if np.linalg.norm(np.array([tgt["SSX"], tgt["SSY"]], dtype=float)) > imageShiftLimit - alignLimit:
                log(f"WARNING: Target [{str(i + 1).zfill(3)}] is too close to the image shift limit. This target will we skipped.")
                skip = True

            if skip: 
                position[-1][0]["skip"] = True
                position[-1].append(copy.deepcopy(position[-1][0]))
                position[-1].append(copy.deepcopy(position[-1][0]))
                skippedTgts += 1
                continue

            if tiltTargets == 0:
                tiltScaling = np.cos(np.radians(pretilt * np.cos(np.radians(rotation)) + startTilt)) / np.cos(np.radians(pretilt * np.cos(np.radians(rotation)))) # stretch shifts from 0 tilt to startTilt
            else:
                tiltScaling = 1
            log(f"DEBUG: Tilt scaling to start tilt [{startTilt}]: {tiltScaling}")

            sem.ImageShiftByMicrons(float(tgt["SSX"]), float(tgt["SSY"]) * tiltScaling)             # apply relative shifts to find out absolute IS after realign to item
            if (previewAli or viewAli):                                                             # adds initial dose, but makes sure start tilt image is on target
                if alignToP:
                    is_x, is_y, *_ = sem.ReportImageShift()
                    if beamTiltComp:
                        saveZeroImageShiftDefocusXLens()
                        doLafis(is_x,is_y)
                        ronchi_before_preview_align(f"target {i + 1} preview alignment (alignToP)", pos=i, pn=0)
                    if beamTiltComp:
                        restoreLafis()
                    x, y, binning, exp, *_ = sem.ImageProperties("P")
                    sem.SetExposure("V", exp)
                    sem.SetBinning("V", int(binning))
                    is_x, is_y, *_ = sem.ReportImageShift()
                    sem.GoToLowDoseArea("V")
                    sem.SetImageShift(0, 0)
                    sem.SetImageShift(is_x, is_y)
                    sem.V()
                    sem.CropCenterToSize("A", int(x), int(y))
                    alignTo("P", debug)
                    sem.RestoreCameraSet("V")
                else:
                    if "viewfile" in tgt.keys() and viewAli and i != 0:                             # skip for tracking target since it was already aligned after tilt to startTilt   
                        sem.ReadOtherFile(0, "O", tgt["viewfile"])                                  # reads view file for first AlignTo instead
                        is_x, is_y, *_ = sem.ReportImageShift()
                        sem.GoToLowDoseArea("V")
                        sem.SetImageShift(0, 0)
                        sem.SetImageShift(is_x, is_y)
                        sem.V()
                        alignTo("O", debug)
                        ASX, ASY = sem.ReportAlignShift()[4:6]
                        log(f"Target alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")    
                    if "tgtfile" in tgt.keys() and previewAli:                
                        sem.ReadOtherFile(0, "O", tgt["tgtfile"])                                   # reads tgt file for first AlignTo instead
                        is_x, is_y, *_ = sem.ReportImageShift()
                        sem.GoToLowDoseArea("R")
                        sem.SetImageShift(0, 0)
                        sem.SetImageShift(is_x, is_y)
                        if beamTiltComp:
                            saveZeroImageShiftDefocusXLens()
                            doLafis(is_x,is_y)
                        ronchi_before_preview_align(f"target {i + 1} preview alignment (tgtfile)", pos=i, pn=0)
                        tem.L()
                        if beamTiltComp:
                            restoreLafis()
                        alignTo("O", debug)
                        AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
                        log(f"Target alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
                    elif "viewfile" in tgt.keys() and previewAli:
                        # Use align between mags to align preview image to view image
                        if not viewAli:
                            #sem.GoToLowDoseArea("V")                                                # If ReadOtherFile while in Record, pixel size of Record is used and AlignBetweenMags fails (seems to be fixed in 4.2beta from 14.08.2024)
                            sem.ReadOtherFile(0, "O", tgt["viewfile"])                              # reads view file for first AlignTo instead
                        else:
                            # If View image was already aligned, take new centered View image at startTilt and use as reference instead
                            is_x, is_y, *_ = sem.ReportImageShift()
                            sem.GoToLowDoseArea("V")
                            sem.SetImageShift(0, 0)
                            sem.SetImageShift(is_x, is_y)
                            sem.V()
                            sem.Copy("A", "O")
                        # Check defocus offset
                        is_x, is_y, *_ = sem.ReportImageShift()
                        sem.GoToLowDoseArea("R")                                                    # Switch to R before applying defocus offset to not mess with potential mP/nP offsets between View and Rec
                        sem.SetImageShift(0, 0)
                        sem.SetImageShift(is_x, is_y)
                        defocus_offset = max(-10, sem.ReportLDDefocusOffset("V"))
                        if defocus_offset != 0:
                            sem.ChangeFocus(defocus_offset)                                             # Higher defocus for better correlation, but max at 10 to avoid major distortions
                        if beamTiltComp:
                            doLafis(is_x,is_y)
                        add_lpp_meta_to_next_mdoc()
                        ronchi_before_preview_align(f"target {i + 1} preview alignment (view to Record)", pos=i, pn=0)
                        sem.L()
                        if beamTiltComp:
                            restoreLafis()
                        sem.AlignBetweenMags("O", -1, -1, -1)
                        AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
                        if defocus_offset != 0:
                            sem.ChangeFocus(-defocus_offset)                                            # Reset focus
                        log(f"Target alignment (Pv2V) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")           

                    # Save preview image as new reference
                    if refFromPreview:
                        sem.OpenNewFile(os.path.splitext(tgt["tgtfile"])[0] + "_tempref.mrc")
                        sem.S()
                        sem.CloseFile()
                        position[-1][0]["ISXali"] = AISX                                            # Save shifts to real reference
                        position[-1][0]["ISYali"] = AISY 

                sem.GoToLowDoseArea("R")
            ISXset, ISYset, *_ = sem.ReportImageShift()
            SSX, SSY = sem.ReportSpecimenShift()
            sem.SetImageShift(ISX0, ISY0)                                                           # reset IS to center position    

            z0_ini = np.tan(np.radians(pretilt)) * (np.cos(np.radians(rotation)) * float(tgt["SSY"]) - np.sin(np.radians(rotation)) * float(tgt["SSX"]))
            correctedFocus = positionFocus - z0_ini * np.cos(np.radians(startTilt)) - float(tgt["SSY"]) * np.sin(np.radians(startTilt))

            position[-1][0]["SSX"] = float(SSX)
            position[-1][0]["SSY"] = float(SSY)
            position[-1][0]["focus"] = correctedFocus
            position[-1][0]["z0"] = z0_ini                                                          # offset from eucentric height (will be refined during collection)
            position[-1][0]["n0"] = float(tgt["SSY"])                                               # offset from tilt axis
            position[-1][0]["ISXset"] = float(ISXset)
            position[-1][0]["ISYset"] = float(ISYset)

            position[-1].append(copy.deepcopy(position[-1][0]))                                     # plus and minus branch start with same values
            position[-1].append(copy.deepcopy(position[-1][0]))

            position[-1][1]["n0"] -= taOffsetPos
            position[-1][2]["n0"] -= taOffsetNeg

            positionFocus += stepDefocus                                                            # adds defocus step between targets and resets to initial defocus if minDefocus is surpassed
            if positionFocus > minFocus0: positionFocus = focus0

    ### Start tilt
        log("Start tilt series...", style=1)

        dewarFillTime = 0
        maxProgress = total_tilt_steps * (len(position) - skippedTgts)
        resumePercent = 0
        startTime = sem.ReportClock()
        lastSlitCheck = startTime

        geo = [[], [], []]
        tiltStepCounter = 0
        schedule_start = 0
        resumePN = 0
        resumePlus = startTilt
        resumeMinus = startTilt
        posResumed = -1

    ### Recovery attempt
    else:
        if realign:
            sem.MoveToNavItem(navID)
            if alignToP:
                x, y, binning, exp, *_ = sem.ImageProperties("P")
                sem.SetExposure("V", exp)
                sem.SetBinning("V", int(binning))
                is_x, is_y, *_ = sem.ReportImageShift()
                sem.GoToLowDoseArea("V")
                sem.SetImageShift(0, 0)
                sem.SetImageShift(is_x, is_y)
                sem.V()
                sem.CropCenterToSize("A", int(x), int(y))
                alignTo("P", debug)
                sem.RestoreCameraSet("V")
            else:
                sem.RealignToOtherItem(navID, 1)
        position = []
        skippedTgts = 0
        for pos in range(len(targets)):
            position.append([{},{},{}])
            for i in range(2):
                position[-1][i+1]["SSX"] = float(savedRun[pos][i]["SSX"])
                position[-1][i+1]["SSY"] = float(savedRun[pos][i]["SSY"])
                position[-1][i+1]["focus"] = float(savedRun[pos][i]["focus"])
                position[-1][i+1]["z0"] = float(savedRun[pos][i]["z0"])
                position[-1][i+1]["n0"] = float(savedRun[pos][i]["n0"])
                if savedRun[pos][i]["shifts"] != "" and not realign:
                    position[-1][i+1]["shifts"] = [float(shift) for shift in savedRun[pos][i]["shifts"].split(",")]
                else:
                    position[-1][i+1]["shifts"] = []
                if savedRun[pos][i]["angles"] != "" and not realign:
                    position[-1][i+1]["angles"] = [float(angle) for angle in savedRun[pos][i]["angles"].split(",")]
                else:
                    position[-1][i+1]["angles"] = []
                position[-1][i+1]["ISXset"] = float(savedRun[pos][i]["ISXset"])
                position[-1][i+1]["ISYset"] = float(savedRun[pos][i]["ISYset"])
                position[-1][i+1]["ISXali"] = float(savedRun[pos][i]["ISXali"])
                position[-1][i+1]["ISYali"] = float(savedRun[pos][i]["ISYali"])
                position[-1][i+1]["dose"] = float(savedRun[pos][i]["dose"])
                position[-1][i+1]["sec"] = int(savedRun[pos][i]["sec"])
                position[-1][i+1]["skip"] = True if savedRun[pos][i]["skip"] == "True" or targets[pos]["skip"] == "True" else False

            sem.AreaForCumulRecordDose(pos + 1)                                                     # set dose accumulator to highest recorded prior dose
            sem.AccumulateRecordDose(max(position[-1][1]["dose"], position[-1][2]["dose"]))

            if targets[pos]["skip"] == "True":
                skippedTgts += 1

        tilt_schedule = build_tilt_schedule(tiltScheme, startTilt, minTilt, maxTilt, step, groupSize)
        total_tilt_steps = len(tilt_schedule)
        log(f"Tilt series: {total_tilt_steps} angles ({tilt_schedule[0]} to {tilt_schedule[-1]} deg)")

        posResumed = resume["pos"] + 1
        tiltStepCounter = resume["sec"]
        schedule_start = resume["sec"]

        resumePlus = startTilt
        resumeMinus = startTilt
        resumePN = 1
        for idx in range(schedule_start):
            t = tilt_schedule[idx]
            if t >= startTilt:
                resumePlus = t
            else:
                resumeMinus = t
        if schedule_start < total_tilt_steps:
            resumePN = 1 if tilt_schedule[schedule_start] >= startTilt else 2
        elif schedule_start > 0:
            resumePN = 1 if tilt_schedule[schedule_start - 1] >= startTilt else 2

        dewarFillTime = 0
        maxProgress = total_tilt_steps * (len(position) - skippedTgts)
        # progress = collected images * (positions - skipped positions) + current position - skipped positions scaled assuming homogeneous distribution of skipped positions
        progress = resume["sec"] * (len(position) - skippedTgts) + resume["pos"] - skippedTgts * resume["pos"] / len(position)
        resumePercent = round(100 * (progress / maxProgress), 1)

        sem.GoToLowDoseArea("R")
        apply_nav_start_defocus(item_index, "recovery")
        origMag, *_ = sem.ReportMag()
        s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
        is2ssMatrix = np.array(sem.ISToSpecimenMatrix(0)).reshape((2, 2))
        ss2isMatrix = np.array(sem.SpecimenToISMatrix(0)).reshape((2, 2))
        camX, camY, *_ = sem.CameraProperties()
        c2ssMatrix = np.array(sem.CameraToSpecimenMatrix(0)).reshape((2, 2))
        ss2cMatrix = np.array(sem.SpecimenToCameraMatrix(0)).reshape((2, 2))
        if debug:
            log("DEBUG: Conversion matrices:")
            log(f"    Stage to Specimen: {s2ssMatrix}", color=1)
            log(f"    IS to Specimen: {is2ssMatrix}", color=1)
            log(f"    Image to Specimen: {c2ssMatrix}", color=1)
            log(f"    Specimen to Camera: {ss2cMatrix}", color=1)

        focus0 = (position[0][1]["focus"] + position[0][2]["focus"]) / 2                            # get estimate for original microscope focus value by taking average of both branches of tracking target

        startTime = sem.ReportClock()
        lastSlitCheck = startTime


    ### Tilt series
    run_tilt_series(schedule_start)

    ### Finish
    sem.ClearStatusLine(0)
    if trackMag > 0:
        sem.RestoreLowDoseParams("R")                                                               # restore record mag before script just in case
    sem.TiltTo(0)
    sem.SetDefocus(focus0)
    sem.SetImageShift(0, 0)
    sem.CloseFile()
    _reset_ronchi_to_session_start("after nav item acquisition")
    updateTargets(runFileName, targets)

    # Format final tilt stacks
    if delFinalStack:
        for target in targets:
            if checkFrames(target["tsfile"]):
                os.remove(target["tsfile"])
                log(f"NOTE: {target['tsfile']} was deleted. Please use saved frames to generate the tilt series.")
    else:
        if sortByTilt or binFinalStack > 1:
            for target in targets:
                if sortByTilt:
                    sortTS(target["tsfile"])
                if binFinalStack > 1:
                    binStack(target["tsfile"], binFinalStack)

    totalTime = round(sem.ReportClock() / 60, 1)
    perTime = round(totalTime / len(position), 1)
    if recoverInput == 1:
        perTime = f"since recovery: {perTime}"
    log(datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
    log(f"##### All tilt series completed in {totalTime} min ({perTime} min per tilt series) #####", color=3, style=1)
    sem.SaveLog()
    return True

def dumpVars(filename):
    output = "# PACEtomo settings from " + datetime.now().strftime("%d.%m.%Y %H:%M:%S") + "\n"
    save = False
    for var in globals():                                                                       # globals() is ordered by creation, start and end points might have to be adjusted if script changes
        if var == "sem":                                                                        # first var after settings vars
            break
        if save:
            output += var + " = " + str(globals()[var]) + "\n"
        if var == "SEMflush":                                                                   # last var before settings vars
            save = True
    with open(filename + "_settings.txt", "w") as f:
        f.write(output)

######## END FUNCTIONS ########

# Adjust user settings
sem.SetProperty("ImageShiftLimit", imageShiftLimit)
tiltLimit = sem.ReportProperty("MaximumTiltAngle")

sem.SetUserSetting("DriftProtection", 1)
sem.SetUserSetting("ShiftToTiltAxis", 1)
sem.SetNewFileType(0)                                                                           # set file type to mrc in case user changed default file type
sem.SetFrameBaseName(0, 1, 0, "PACEtomo_setup")                                                 # change frame name at start to avoid overwriting in case sets other than Record save frames

# Warnings
log(f"DEBUG: Tilt limit is: {tiltLimit}")
if (maxTilt > tiltLimit or minTilt < -tiltLimit) and sem.IsVariableDefined("warningTiltAngle") == 0:
    sem.Pause("WARNING: Tilt angles go beyond +/- 70 degrees. Most stage limitations do not allow for symmetrical tilt series with these values!")
    sem.SetPersistentVar("warningTiltAngle", "")

if int(sem.ReportAxisPosition("F")[0]) != 0 and sem.IsVariableDefined("warningFocusArea") == 0:
    sem.Pause("WARNING: Position of Focus area is not 0! Please set it to 0 to autofocus on the tracking target!")
    sem.SetPersistentVar("warningFocusArea", "")

tiltAxisOffset = sem.ReportTiltAxisOffset()[0]
if float(tiltAxisOffset) == 0 and sem.IsVariableDefined("warningTAOffset") == 0:
    sem.Pause("WARNING: No tilt axis offset was set! Please run the PACEtomo_measureOffset script to determine appropiate tilt axis offset.")
    sem.SetPersistentVar("warningTAOffset", "")

# Saving SerialEM setup
sem.SaveSettings()
sem.SaveNavigator()

sem.SuppressReports()
if beamTiltComp:                                                                                # check if there is a calibration saved, throws error if not
    sem.ReportComaVsISmatrix()
if tgtPattern:                                                                                  # initialize in case tgts file contains values
    vecA0 = vecA1 = vecB0 = vecB1 = size = None



### Multi-nav batch entry
log(f"PACEtomo Version {versionPACE}", color=5, style=1)
sem.ProgramTimeStamps()

items_to_run = resolve_nav_items()
validate_nav_lists()
checkRonchigramSetup()

start_idx = batch_resume_start_index(items_to_run)
if start_idx >= len(items_to_run):
    log("NOTE: All nav items in list appear complete.")
    sem.Exit()

batch_recover_accepted = False
if start_idx < len(items_to_run):
    curDir_batch = sem.ReportDirectory()
    tf_path = get_nav_tgts_path(items_to_run[start_idx], curDir_batch)
    if tf_path and nav_run_status(tf_path) == "incomplete":
        nav_idx_resume = items_to_run[start_idx]
        batch_recover_accepted = sem.YesNoBox(
            f"The target file for nav item {nav_idx_resume} contains recovery data. "
            "Do you want to attempt to continue the acquisition? Tracking accuracy might be impacted."
        ) == 1
        if not batch_recover_accepted:
            sem.AllowFileOverwrite(1)

for item_index in range(start_idx, len(items_to_run)):
    nav_idx = items_to_run[item_index]
    log(f"===== Nav item {nav_idx} ({item_index + 1}/{len(items_to_run)}) =====", style=1)
    ok = run_one_nav_item(
        nav_idx,
        item_index,
        batch_recover=(item_index == start_idx),
        batch_recover_accepted=batch_recover_accepted,
    )
    if not ok:
        log(f"WARNING: Nav item {nav_idx} was skipped due to an error.")
    if item_index < len(items_to_run) - 1:
        reset_scope_between_nav_items()

sem.Exit()
