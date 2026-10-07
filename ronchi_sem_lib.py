#!Python
# need sem
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
import ronchi_lib
import numpy as np

debug = False

hasXLens = True
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
ronchiMeasureCount = 0
########## END Ronchigram settings ##########

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

def add_lpp_meta_to_next_mdoc():
    for k,v in (
            ('ImageDistanceOffset', sem.ReportImageDistanceOffset()),
        ):
        v_str = '%.12f' % (float(v))
        sem.AddToNextFrameStackMdoc(k, v_str)

def checkRonchigramSetup():
    global ronchiStartXLensX, ronchiStartXLensY, ronchiStartC3Offset, doRonchigram, ronchiC3Offset
    if ronchiStartXLensX is None:
        ronchiStartXLensX, ronchiStartXLensY = [float(v) for v in sem.ReportXLensDeflector(2)[:2]]
    if ronchiStartC3Offset is None:
        ronchiStartC3Offset = float(sem.ReportImageDistanceOffset())
    if not ronchiC3Offset:
        sem.Pause('Please set C3 offset to where you can clearly see the global xLPP center')
    
        ronchiC3Offset = float(sem.ReportImageDistanceOffset()) - ronchiStartC3Offset

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
        log('WARNING: illuminated area not available', e)
        illuminated_area = None
    if ronchi_c3_offset is not None:
        ronchi_c3_offset = float(ronchi_c3_offset)
    else:
        ronchi_c3_offset = None
    kv_pairs = {'mag': int(mag), 'spot_size': int(spot_size), 'illuminated_area': illuminated_area,
            'cam_binning': cam_binning, 'ronchi_binning': ronchi_binning, 'ronchi_c3_offset':ronchi_c3_offset}
    return kv_pairs

###############
# Acquiring Ronchigram
##############
def acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=''):
    sem.SetImageDistanceOffset(trial_offset_baseline + ronchi_offset)
    # We don't care about saving now.
    #saved_basename = _set_ronchi_trial_frame_basename()
    my_image = None
    try:
        sem.Delay(ronchiDelay, "s")
        add_lpp_meta_to_next_mdoc()
        # acquire with preset parameters
        getattr(sem,sem_acquire_preset)()
        my_image = np.asarray(sem.bufferImage("A"))
    finally:
        sem.SetImageDistanceOffset(trial_offset_baseline)
        #_restore_frame_basename(saved_basename)
    if pass_label:
        log(f"Ronchigram{pass_label}: Trial image acquired.")
    return my_image

def _reset_ronchi_xlens_if_out_of_tolerance(lens_index=2):
    pass

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

def _acquire_ronchi_trial(trial_offset_baseline, pass_label=""):
    """Trial ronchigram at C3 imaging offset; restore baseline offset after shot."""
    _reset_ronchi_xlens_if_out_of_tolerance()
    is_x, is_y, *_ = sem.ReportImageShift()
    sem.GoToLowDoseArea("T")
    sem.SetImageShift(0, 0)
    sem.SetImageShift(is_x, is_y)
    #saved_basename = _set_ronchi_trial_frame_basename()
    my_image = acquire_ronchi_image(trial_offset_baseline, ronchiC3Offset, 'T',pass_label)
    #_restore_frame_basename(saved_basename)
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
    #global ronchiMeasureCount
    image = _acquire_ronchi_trial(c3_baseline_offset, pass_label=pass_label)
    result = _analyze_ronchi_image(image)
    #ronchiMeasureCount += 1
    _log_ronchi_ks(result, pass_label=pass_label)
    _log_ronchi_phases(result, pass_label=pass_label)
    return result

# Corrections
def applyRonchigramXtiltCorrection(correction_x, correction_y, lens_index=2):
    """Apply calculated shifts to the X lens deflector."""
    if not hasXLens:
        return
    xtX, xtY = sem.ReportXLensDeflector(lens_index)
    sem.SetXLensDeflector(lens_index, xtX + correction_x, xtY + correction_y)

def _apply_ronchi_phase(result, pass_label=""):
    _log_ronchi_phases(result, pass_label=pass_label)
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
    if abs(result["c3_correction"]) <= min_err:
        return False
    new_offset = applyRonchigramC3Correction(result["c3_correction"], c3_baseline_offset)
    return True

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
