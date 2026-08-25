#!Python
# need sem
import serialem as sem
import ronchi_lib
import numpy as np

hasXLens = True
doRonchigram       = True
ronchiBaseSuffix   = "_ronchi"         # appended to active frame base name for Trial saves only, then restored
ronchiC3Offset     = -20          # added to ReportImageDistanceOffset before Trial shot
ronchiDelay        = 1.0          # seconds after C3 offset change
ronchiBinning      = 2
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

def acquire_ronchi_image(trial_offset_baseline, ronchi_offset, sem_acquire_preset='T',pass_label=''):
    sem.SetImageDistanceOffset(trial_offset_baseline + ronchi_offset)
    # We don't care about saving now.
    #saved_basename = _set_ronchi_trial_frame_basename()
    try:
        sem.Delay(ronchiDelay, "s")
        add_lpp_meta_to_next_mdoc()
        # acquire with preset parameters
        getattr(sem,sem_acquire_preset)()
    finally:
        sem.SetImageDistanceOffset(trial_offset_baseline)
        #_restore_frame_basename(saved_basename)
    if pass_label:
        log(f"Ronchigram{pass_label}: Trial image acquired.")
    return np.asarray(sem.bufferImage("A"))
