#!Python
##############################################################################
# Ronchigram - microscope aligner (SerialEM: Trial acquire, analyze, apply)
#
# Shared by PACEtomo.py, PACEtomo_selectTargets.py and the calibration scripts.
# Image analysis (numpy only) is in ronchi_lib.py.
#
# Settings below are defaults for standalone use (ronchi_cal.py, lafis_cal.py,
# ronchi_test.py ...).  PACEtomo scripts keep their own settings in their
# header (so _bset overrides and the settings dump keep working) and push
# them here with configure(globals(), log_fn=log) before use.
#
# Public API: configure(), getOpticsKey(), readCalibrations(),
#             checkRonchigramSetup(), acquire_ronchi_image(),
#             doRonchiCorrection(), doRonchigramCorrection(),
#             reset_ronchi_measure_count(), add_lpp_meta_to_next_mdoc()
##############################################################################
import sys
import platform
if platform.system() == 'Windows':
    is_simu = False
    required_sys_paths = [    'C:\Program Files\SerialEM\PythonModules',
                    'C:\\Users\\VALUEDGATANCUSTOMER\\Desktop\\anchi\\PACEtomo_LPP',
                ]
    for p in required_sys_paths:
        if p not in sys.path:
            sys.path.insert(0,p)
    import serialem as sem
else:
    is_simu = True
    print('testing on Mac/Linux with simulator')
    import sem_simulator as sem

import numpy as np
import ronchi_lib
import cal_util

########## Ronchigram scope constants (not taken from callers by configure) ##########
ronchiBinning      = 32
ronchiPixelSize    = 0.98e-4 * 2 # um (unbinned; multiplied by binning in analysis)
ronchiCorrMatrix   = [[0.212, 1.28], [1.22, -0.243]]  # phase-to-deflector coupling, scaled by 1e-5

########## Ronchigram settings (defaults; overridden by configure) ##########
debug = False

hasXLens = True
doRonchigram       = True
ronchiBaseSuffix   = "_ronchi"         # appended to active frame base name for Trial saves only, then restored
ronchiC3Offset     = -20          # added to ReportImageDistanceOffset before Trial shot
ronchiDelay        = 1.0          # seconds after C3 offset change
ronchiTargetPhaseA = -1.93941993           # vertical laser (rad)
ronchiTargetPhaseB = 1.67658165        # horizontal laser (rad)
ronchiCorrectKs    = [[9.303, -0.662] ,  [0.856 ,8.680]]
ronchiPeakRadius   = 100
ronchiCorrectC3    = True         # apply C3 correction from mean ks error (diagonal fringe spacing)
ronchiC3CorrectionFactor = 20 / 9.1  # um offset per um^-1 mean ks error
ronchiCorrectConstantIterNumber  = 0         # measure and correct c3 and phase on a constant number of iteration.  ronchiMinErrForC3Correction, ronchiMinErrForC3CorrectionRedo, and redo_ronchi_after_C3 are ignored when this is set to a number larger than 0.
ronchiMinErrForC3Correction   = 0.3          # apply C3 on 1st Trial only if |c3 correction| exceeds this (um)
ronchiMinErrForC3CorrectionRedo = 0.5        # apply C3 on 2nd Trial only if |c3 correction| exceeds this (um)
redo_ronchi_after_C3 = True       # up to 3 Trials: 1st C3, 2nd optional C3 + 3rd phase-only if 2nd C3 applied
ronchiXLensTolerance = 0.000125     # reset XLensDeflector(2) to start if |x-x0| or |y-y0| exceeds this
ronchiXLensOffset  = (0.0, 0.0)     # XLens correction made by others (e.g. PACEtomo LAFIS); excluded from tolerance check
########## END Ronchigram settings ##########

# Names configure() accepts from a caller namespace
SETTING_NAMES = (
    "debug", "hasXLens", "doRonchigram", "ronchiBaseSuffix", "ronchiC3Offset", "ronchiDelay",
    "ronchiTargetPhaseA", "ronchiTargetPhaseB",
    "ronchiCorrectKs", "ronchiPeakRadius", "ronchiCorrectC3",
    "ronchiC3CorrectionFactor", "ronchiCorrectConstantIterNumber",
    "ronchiMinErrForC3Correction", "ronchiMinErrForC3CorrectionRedo",
    "redo_ronchi_after_C3", "ronchiXLensTolerance", "ronchiXLensOffset",
)

########## Session state (owned by this module) ##########
ronchiStartXLensX = None          # set from ReportXLensDeflector(2) in checkRonchigramSetup
ronchiStartXLensY = None
ronchiStartC3Offset = None        # set from ReportImageDistanceOffset in checkRonchigramSetup
ronchiMeasureCount = 0            # incremented in _ronchi_trial_and_analyze to log result in frame stack mdoc


def _default_log(text, color=0, style=0):
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
    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)

_log_fn = _default_log

def log(text, color=0, style=0):
    return _log_fn(text, color, style)


def configure(namespace=None, log_fn=None, **settings):
    """
    Copy caller settings into this module.
    namespace: dict such as the caller's globals(); only names in SETTING_NAMES are taken.
    log_fn: caller's log(text, color, style) so messages go through its formatting.
    settings: explicit overrides, applied after namespace.
    """
    global _log_fn
    g = globals()
    if namespace is not None:
        for k in SETTING_NAMES:
            if k in namespace:
                g[k] = namespace[k]
    for k, v in settings.items():
        if k not in SETTING_NAMES:
            raise KeyError(f"unknown ronchigram setting {k}")
        g[k] = v
    if log_fn is not None:
        _log_fn = log_fn


###############
# Optics and calibrations
##############
def getOpticsKey(ronchi_binning=32, ronchi_c3_offset=None):
    """
    Optics condition that calibrations such as pixel_xt_matrix depend on.
    ronchi_binning: the binning of the ronchi_image before fft is calculated for
    peak finding
    ronchi_c3_offset: ronchi_sem_lib.ronchiC3Offset used to acquire
    the ronchigram.
    This is in ronchi_sem_lib.py instead of cal_util.py to avoid need to
    import serialem in cal_util.py
    """
    mag,*_ = sem.ReportMag()
    spot_size = sem.ReportSpotSize()
    cam_binning = sem.ReportBinning('T')
    try:
        illuminated_area = float(sem.ReportIlluminatedArea())
    except Exception as e:
        # only available on some Thermo Scientific microscopes
        log(f'WARNING: illuminated area not available {e}')
        illuminated_area = None
    if ronchi_c3_offset is not None:
        ronchi_c3_offset = float(ronchi_c3_offset)
    kv_pairs = {'mag': int(mag), 'spot_size': int(spot_size), 'illuminated_area': illuminated_area,
            'cam_binning': cam_binning, 'ronchi_binning': ronchi_binning, 'ronchi_c3_offset':ronchi_c3_offset}
    return kv_pairs

def readCalibrations(sem_mode='T', ronchi_binning=32, ronchi_c3_offset=None):
    """
    Read ronchi_ref_ks and ronchi_ref_phase calibrations at the current optics.
    Look for ronchi_c3_offset first. If none is found, fall back to the newest
    calibration at any C3 offset and warn, since ks depend on the C3 offset.
    Returns (target_phase_a, target_phase_b, correct_ks, calibration_c3_offset).
    """
    sem.GoToLowDoseArea(sem_mode)
    working_dir = sem.ReportDirectory()
    cal_dir, session_name = cal_util.getCalibrationsDir(working_dir)
    optics = getOpticsKey(ronchi_binning=ronchi_binning, ronchi_c3_offset=ronchi_c3_offset)
    ref_correct_ks, cal_c3_offset = cal_util.readCalibration('ronchi_ref_ks', cal_dir, optics)
    if ref_correct_ks is None and ronchi_c3_offset is not None:
        optics = getOpticsKey(ronchi_binning=ronchi_binning, ronchi_c3_offset=None)
        ref_correct_ks, cal_c3_offset = cal_util.readCalibration('ronchi_ref_ks', cal_dir, optics)
        if ref_correct_ks is not None:
            log(f"WARNING: no ronchi_ref_ks calibration at C3 offset {ronchi_c3_offset}; "
                f"using the one at C3 offset {cal_c3_offset}.")
    if ref_correct_ks is None:
        raise RuntimeError(f"No ronchi_ref_ks calibration in {cal_dir} for optics {optics}")
    # phases must come from the same C3 offset as ks
    optics = getOpticsKey(ronchi_binning=ronchi_binning, ronchi_c3_offset=cal_c3_offset)
    ref_phases, cal_c3_offset = cal_util.readCalibration('ronchi_ref_phase', cal_dir, optics)
    if ref_phases is None:
        raise RuntimeError(f"No ronchi_ref_phase calibration in {cal_dir} for optics {optics}")
    return ref_phases[0], ref_phases[1], ref_correct_ks, cal_c3_offset


###############
# Setup
##############
def checkRonchigramSetup():
    """Startup checks when doRonchigram is enabled. Records session start XLens and C3."""
    global ronchiStartXLensX, ronchiStartXLensY, ronchiStartC3Offset
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
    if ronchiCorrectC3:
        if ronchiCorrectConstantIterNumber == 0:
            c3_spec_txt = f"only if |C3 correction| > {ronchiMinErrForC3Correction} um"
        else:
            c3_spec_txt = f"{ronchiCorrectConstantIterNumber:d} iterations, only if |C3 correction| > {ronchiMinErrForC3Correction} um"
        log(f"NOTE: Ronchigram C3 correction enabled (factor {ronchiC3CorrectionFactor:.4f}, {c3_spec_txt}).")
    if ronchiCorrectConstantIterNumber > 0:
        log(f"NOTE: Ronchigram uses {ronchiCorrectConstantIterNumber:d} Trials, each applying C3 and phase corrections.")
    elif redo_ronchi_after_C3:
        log(
            "NOTE: Ronchigram may use up to 3 Trials before Record: "
            "1st C3 (if above threshold), 2nd phase or C3, 3rd phase-only if 2nd C3 applied."
        )
        log(
            f"NOTE: 2nd-Trial C3 threshold |correction| > {ronchiMinErrForC3CorrectionRedo} um "
            f"(1st-Trial threshold {ronchiMinErrForC3Correction} um)."
        )


###############
# Frame stack mdoc
##############
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

def reset_ronchi_measure_count():
    global ronchiMeasureCount
    for c0 in range(ronchiMeasureCount):
        # Enter default value for each key so that
        # it does not retain old value in the next mdoc.
        # SerialEM's default behavior is such that key-value pair
        # persists until the program is closed.
        c1 = c0+1 # base 1
        sem.AddToNextFrameStackMdoc(f'RonchiKs{c1:02d}', '')
        sem.AddToNextFrameStackMdoc(f'RonchiPhases{c1:02d}','')
    # reset value until the next target
    ronchiMeasureCount = 0


###############
# Frame base name for Trial saves
##############
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


###############
# Acquiring Ronchigram
##############
def acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=''):
    """Acquire at baseline + ronchi_offset C3; restore baseline offset after shot. No frame base name change."""
    sem.SetImageDistanceOffset(trial_offset_baseline + ronchi_offset)
    my_image = None
    try:
        sem.Delay(ronchiDelay, "s")
        add_lpp_meta_to_next_mdoc()
        # acquire with preset parameters
        getattr(sem,sem_acquire_preset)()
        my_image = np.asarray(sem.bufferImage("A"))
    finally:
        sem.SetImageDistanceOffset(trial_offset_baseline)
    if pass_label:
        log(f"Ronchigram{pass_label}: Trial image acquired.")
    return my_image

def _reset_ronchi_xlens_if_out_of_tolerance(lens_index=2):
    """Reset XLensDeflector to session start if phase corrections have walked it too far."""
    if not hasXLens or ronchiStartXLensX is None or ronchiStartXLensY is None:
        return
    offset_x, offset_y = ronchiXLensOffset
    xtX, xtY = sem.ReportXLensDeflector(lens_index)
    xtX, xtY = float(xtX), float(xtY)
    xtX_delta = xtX - offset_x - ronchiStartXLensX
    xtY_delta = xtY - offset_y - ronchiStartXLensY
    if (abs(xtX_delta) > ronchiXLensTolerance
            or abs(xtY_delta) > ronchiXLensTolerance):
        log(
            f"WARNING: Ronchigram X lens deflector ({xtX:.6f}, {xtY:.6f}) beyond tolerance "
            f"{ronchiXLensTolerance} from start ({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}); resetting."
        )
        sem.SetXLensDeflector(lens_index, ronchiStartXLensX + offset_x, ronchiStartXLensY + offset_y)

def _acquire_ronchi_trial(trial_offset_baseline, pass_label=""):
    """Trial ronchigram at C3 imaging offset; restore baseline offset and frame base name after shot."""
    _reset_ronchi_xlens_if_out_of_tolerance()
    is_x, is_y, *_ = sem.ReportImageShift()
    sem.GoToLowDoseArea("T")
    sem.SetImageShift(0, 0)
    sem.SetImageShift(is_x, is_y)
    saved_basename = _set_ronchi_trial_frame_basename()
    try:
        my_image = acquire_ronchi_image(trial_offset_baseline, ronchiC3Offset, 'T', pass_label)
    finally:
        _restore_frame_basename(saved_basename)
    return np.asarray(my_image)

###############
# Analyzing Ronchigram
##############
def _analyze_ronchi_image(image):
    return ronchi_lib.analyze_ronchigram(
        image, ronchiPixelSize, ronchiBinning, ronchiTargetPhaseA, ronchiTargetPhaseB,
        ronchiCorrectKs, peak_radius=ronchiPeakRadius, corr_matrix=ronchiCorrMatrix,
        c3_correction_factor=ronchiC3CorrectionFactor,
    )


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
    sem.AddToNextFrameStackMdoc(f'RonchiKs{ronchiMeasureCount:02d}', f"{np.array2string(result['ks'].ravel(), precision=4)[1:-1]}")

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
    sem.AddToNextFrameStackMdoc(f'RonchiPhases{ronchiMeasureCount:02d}',f"{phases[0]:.3f}    {phases[1]:.3f}")

def _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=""):
    """Acquire Trial ronchigram and analyze. Returns analysis result dict."""
    global ronchiMeasureCount
    image = _acquire_ronchi_trial(c3_baseline_offset, pass_label=pass_label)
    result = _analyze_ronchi_image(image)
    ronchiMeasureCount += 1
    _log_ronchi_ks(result, pass_label=pass_label)
    _log_ronchi_phases(result, pass_label=pass_label)
    return result

###############
# Corrections
##############
def applyRonchigramXtiltCorrection(correction_x, correction_y, lens_index=2):
    """Apply calculated shifts to the X lens deflector."""
    if not hasXLens:
        return
    xtX, xtY = sem.ReportXLensDeflector(lens_index)
    sem.SetXLensDeflector(lens_index, xtX + correction_x, xtY + correction_y)

def _apply_ronchi_phase(result, pass_label=""):
    log(f"Ronchigram{pass_label}: applying deflector dX={result['correction_x']:.3e} dY={result['correction_y']:.3e}")
    applyRonchigramXtiltCorrection(result["correction_x"], result["correction_y"])

def applyRonchigramC3Correction(c3_correction, baseline_offset):
    """Apply C3 correction to ImageDistanceOffset relative to pre-ronchigram offset."""
    new_offset = baseline_offset + c3_correction
    sem.SetImageDistanceOffset(new_offset)
    return new_offset

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

def doRonchiCorrectionWithVariableIters(c3_baseline_offset):
    """1st Trial C3; optional 2nd Trial C3 or phase; 3rd Trial phase-only if 2nd C3 applied."""
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
    return result

def doRonchiCorrectionWithConstantIters(c3_baseline_offset, number_of_iters=3):
    """Apply C3 (if above threshold) and phase correction on every one of number_of_iters Trials."""
    result = None
    for itr in range(number_of_iters):
        pass_label = f" ({itr + 1}/{number_of_iters})"
        result = _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=pass_label)
        _try_apply_ronchi_c3(result, c3_baseline_offset, pass_label=pass_label)
        _apply_ronchi_phase(result, pass_label=pass_label)
        # next correction is relative to the C3 just applied
        c3_baseline_offset = float(sem.ReportImageDistanceOffset())
    return result

def doRonchiCorrection(c3_baseline_offset):
    """
    Trial + analyze + C3 and/or phase correction starting from c3_baseline_offset.
    On failure, log, reset C3 to session start and re-raise to the caller.
    Returns the last analysis result dict.
    """
    try:
        if ronchiCorrectConstantIterNumber == 0:
            result = doRonchiCorrectionWithVariableIters(c3_baseline_offset)
        elif ronchiCorrectConstantIterNumber > 0:
            result = doRonchiCorrectionWithConstantIters(c3_baseline_offset, ronchiCorrectConstantIterNumber)
        else:
            raise ValueError(f'Invalid ronchiCorrectConstantIterNumber={ronchiCorrectConstantIterNumber:d}')
        if debug:
            log(f"DEBUG: Ronchigram pixel size {ronchiPixelSize} um, peak radius {ronchiPeakRadius} px")
        return result
    except Exception as e:
        log(f"WARNING: Ronchigram analysis failed: {e}. Continuing without laser correction.")
        if ronchiStartC3Offset is not None:
            sem.SetImageDistanceOffset(ronchiStartC3Offset)
        raise

def doRonchigramCorrection(c3_baseline_offset=None):
    """
    Full Trial correction at the current position, then return to the Record area.
    c3_baseline_offset: None uses the current ImageDistanceOffset.
    Errors are logged and swallowed. Returns the last analysis result dict or None.
    """
    if not doRonchigram:
        return None
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
    if c3_baseline_offset is None:
        c3_baseline_offset = float(sem.ReportImageDistanceOffset())
    result = None
    try:
        result = doRonchiCorrection(c3_baseline_offset)
    except Exception:
        pass  # already logged in doRonchiCorrection
    sem.GoToLowDoseArea("R")
    return result
