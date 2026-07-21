#!Python
# ===================================================================
# ScriptName     test_beam_tilt_defocus
# Purpose:       Autofocus to a single target defocus using the selected
#                measure method, then compare with CtfFind at CTF X-tilt.
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import serialem as sem
import PACEtomo_beamTiltDefocus as btdef

############ SETTINGS ############

# ctf | beam_tilt (physics+calibration) | beam_tilt_sem (sem.G(-1)+calibration)
defocusMethod = "beam_tilt"
# False on scopes without XLensDeflector: skip all XLens Report/Set/Restore.
hasXLens = True

tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
defocus_tilt_correction = beam_tilt_correction
beam_tilt_axis = "x"
spherical_aberration_mm = 2.7

target_defocus_um = -4.0
autofocus_tolerance_um = 0.05
max_autofocus_cycles = 3

ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2
ctf_max_attempts = 3

xtilt_lens_index = 2
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867
beam_tilt_xtilt_x = 0.0
beam_tilt_xtilt_y = 0.0

# Scaling JSON for beam_tilt or beam_tilt_sem; leave "" for uncorrected base.
calibration_file = r""

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


btdef.configure(sem_module=sem, logger=echo, has_x_lens=hasXLens)
if calibration_file:
    btdef.calibration_file = calibration_file
    btdef._calibration_cache = None


def set_xtilt(x, y):
    if not hasXLens:
        return
    sem.SetXLensDeflector(xtilt_lens_index, float(x), float(y))


def run_ctffind():
    sem.NoMessageBoxOnError(1)
    try:
        cfind = sem.CtfFind("A", ctf_defocus_lo, ctf_defocus_hi)
    finally:
        sem.NoMessageBoxOnError(0)
    if len(cfind) == 0:
        sem.Echo("ERROR: CtfFind failed.")
        sem.Exit()
    return float(cfind[0]), float(cfind[-1])


def acquire_ctf_reference():
    if not hasXLens:
        sem.L()
        return run_ctffind()
    xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    try:
        set_xtilt(ctf_xtilt_x, ctf_xtilt_y)
        sem.L()
        return run_ctffind()
    finally:
        set_xtilt(float(xtilt[0]), float(xtilt[1]))


def measure_ctf_defocus():
    ctf_defocus, ctf_res = acquire_ctf_reference()
    raw = {
        "measurement_method": "ctf",
        "legacy_defocus_um": ctf_defocus,
        "defocus_um": ctf_defocus,
        "calibration_correction_um": 0.0,
        "ctf_resolution_A": ctf_res,
        "drift_speed_x_nm_per_s": 0.0,
        "drift_speed_y_nm_per_s": 0.0,
    }
    return ctf_defocus, raw


def measure_defocus_with_diagnostics():
    xt_x = beam_tilt_xtilt_x if hasXLens else None
    xt_y = beam_tilt_xtilt_y if hasXLens else None
    if defocusMethod == "beam_tilt":
        return btdef.measure_defocus_with_diagnostics(
            tilt_angle_mrad=tilt_angle_mrad,
            beam_tilt_correction=beam_tilt_correction,
            defocus_tilt_correction=defocus_tilt_correction,
            xtilt_x=xt_x,
            xtilt_y=xt_y,
            lens_index=xtilt_lens_index,
            beam_tilt_axis=beam_tilt_axis,
            cs_mm=spherical_aberration_mm,
        )
    if defocusMethod == "beam_tilt_sem":
        return btdef.measure_serialEM_defocus_with_diagnostics(
            xtilt_x=xt_x,
            xtilt_y=xt_y,
            lens_index=xtilt_lens_index,
        )
    if defocusMethod == "ctf":
        return measure_ctf_defocus()
    echo(
        f"ERROR: Unknown defocusMethod '{defocusMethod}'. "
        "Use 'ctf', 'beam_tilt', or 'beam_tilt_sem'."
    )
    sem.Exit()


def autofocus_drift_suffix(raw):
    if defocusMethod == "beam_tilt":
        return (
            f"drift=({raw['drift_speed_x_nm_per_s']:.2f}, "
            f"{raw['drift_speed_y_nm_per_s']:.2f}) nm/s"
        )
    return ""


def autofocus_to_target(target):
    """Iterate measure + ChangeFocus until within tolerance."""
    defocus = float("nan")
    raw = {}
    for cycle in range(1, max_autofocus_cycles + 1):
        defocus, raw = measure_defocus_with_diagnostics()
        error = float(target) - defocus
        drift_text = autofocus_drift_suffix(raw)
        suffix = f", {drift_text}" if drift_text else ""
        echo(
            f"Autofocus {cycle}/{max_autofocus_cycles}: "
            f"measured={defocus:.4f} um, target={float(target):.3f} um, "
            f"error={error:.3f} um{suffix}"
        )
        if abs(error) <= autofocus_tolerance_um:
            return defocus, raw
        sem.ChangeFocus(error)
    echo("WARNING: Autofocus did not reach tolerance.")
    return defocus, raw


def echo_final_results(beam_defocus, raw, ctf_defocus, ctf_res):
    base_um = float(
        raw.get("serialEM_defocus_um", raw.get("legacy_defocus_um", beam_defocus))
    )
    correction = float(raw.get("calibration_correction_um", 0.0))
    error_vs_target = beam_defocus - target_defocus_um
    delta_vs_ctf = beam_defocus - ctf_defocus

    echo("================================================")
    echo("FINAL RESULTS")
    echo(f"  defocus method:         {defocusMethod}")
    echo(f"  target defocus:         {target_defocus_um:.4f} um")
    echo(f"  measured (calibrated):  {beam_defocus:.4f} um")
    echo(f"  measured (base):        {base_um:.4f} um")
    echo(f"  calibration correction: {correction:.4f} um")
    if defocusMethod == "beam_tilt" and "cs_term_um" in raw:
        echo(f"  Cs term:                {raw['cs_term_um']:.4f} um")
    if defocusMethod == "beam_tilt_sem":
        echo(f"  measurement:            sem.G(-1)")
    echo(f"  CTF reference:          {ctf_defocus:.4f} um ({ctf_res:.1f} A)")
    echo(f"  error vs target:        {error_vs_target:.4f} um")
    echo(f"  delta calibrated-CTF:   {delta_vs_ctf:.4f} um")
    echo(f"  delta base-CTF:         {base_um - ctf_defocus:.4f} um")
    echo("================================================")


def main():
    sem.SuppressReports()
    echo("##### test_beam_tilt_defocus #####")
    echo(f"defocusMethod={defocusMethod}")
    echo(f"hasXLens={hasXLens}")
    echo(f"calibration_file={btdef.calibration_file or '(none)'}")
    if defocusMethod == "beam_tilt":
        echo(
            f"beam_tilt_correction={beam_tilt_correction}, "
            f"defocus_tilt_correction={defocus_tilt_correction}, "
            f"tilt_angle_mrad={tilt_angle_mrad}"
        )
    echo(f"target defocus={target_defocus_um:.3f} um")
    if hasXLens:
        echo(f"CTF X-tilt=({ctf_xtilt_x:.6f}, {ctf_xtilt_y:.6f})")
        echo(f"measure X-tilt=({beam_tilt_xtilt_x:.6f}, {beam_tilt_xtilt_y:.6f})")
    else:
        echo("NOTE: hasXLens=False; skipping all XLens Report/Set/Restore")

    original_xtilt = None
    if hasXLens:
        original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
        echo(
            f"Saved XLensDeflector({xtilt_lens_index}): "
            f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
        )

    try:
        sem.GoToLowDoseArea("R")
        sem.SetImageShift(0, 0)

        echo("------------------------------------------------")
        echo(f"Autofocus ({defocusMethod})")
        beam_defocus, raw = autofocus_to_target(target_defocus_um)

        echo("------------------------------------------------")
        echo("CTF reference measurement")
        ctf_defocus, ctf_res = acquire_ctf_reference()

        echo_final_results(beam_defocus, raw, ctf_defocus, ctf_res)

    finally:
        if hasXLens and original_xtilt is not None:
            set_xtilt(float(original_xtilt[0]), float(original_xtilt[1]))
            echo(
                f"Restored XLensDeflector({xtilt_lens_index}) to "
                f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
            )

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
