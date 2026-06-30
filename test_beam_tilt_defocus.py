#!Python
# ===================================================================
# ScriptName     test_beam_tilt_defocus
# Purpose:       Beam-tilt autofocus to a single target defocus, then
#                compare the converged value with CtfFind at CTF X-tilt.
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import serialem as sem
import PACEtomo_beamTiltDefocus as btdef

############ SETTINGS (match calibrate_beam_tilt_scaling.py) ############

tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
defocus_tilt_correction = beam_tilt_correction
beam_tilt_axis = "x"

target_defocus_um = -4.0
autofocus_tolerance_um = 0.05
max_autofocus_cycles = 3

ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

xtilt_lens_index = 2
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867
beam_tilt_xtilt_x = 0.0
beam_tilt_xtilt_y = 0.0

# Point at scaling JSON to test calibration; leave "" for physics-only.
calibration_file = r"X:\k3f_leginonframes\p26jun29a\xtilt_calib_test\beam_tilt_scaling_calibration.json"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


btdef.configure(sem_module=sem, logger=echo)
if calibration_file:
    btdef.calibration_file = calibration_file
    btdef._calibration_cache = None


def set_xtilt(x, y):
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
    xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    try:
        set_xtilt(ctf_xtilt_x, ctf_xtilt_y)
        sem.L()
        return run_ctffind()
    finally:
        set_xtilt(float(xtilt[0]), float(xtilt[1]))


def beam_tilt_autofocus(target):
    """Iterate beam-tilt measure + ChangeFocus until within tolerance."""
    defocus = float("nan")
    for cycle in range(1, max_autofocus_cycles + 1):
        defocus, raw = btdef.measure_defocus_with_diagnostics(
            tilt_angle_mrad=tilt_angle_mrad,
            beam_tilt_correction=beam_tilt_correction,
            defocus_tilt_correction=defocus_tilt_correction,
            xtilt_x=beam_tilt_xtilt_x,
            xtilt_y=beam_tilt_xtilt_y,
            lens_index=xtilt_lens_index,
            beam_tilt_axis=beam_tilt_axis,
        )
        error = float(target) - defocus
        echo(
            f"Autofocus {cycle}/{max_autofocus_cycles}: "
            f"measured={defocus:.4f} um, target={float(target):.3f} um, "
            f"error={error:.3f} um, "
            f"drift=({raw['drift_speed_x_nm_per_s']:.2f}, "
            f"{raw['drift_speed_y_nm_per_s']:.2f}) nm/s"
        )
        if abs(error) <= autofocus_tolerance_um:
            return defocus, raw
        sem.ChangeFocus(error)
    echo("WARNING: Beam-tilt autofocus did not reach tolerance.")
    return defocus, raw


def main():
    sem.SuppressReports()
    echo("##### test_beam_tilt_defocus #####")
    echo(f"calibration_file={btdef.calibration_file or '(none)'}")
    echo(
        f"beam_tilt_correction={beam_tilt_correction}, "
        f"defocus_tilt_correction={defocus_tilt_correction}, "
        f"tilt_angle_mrad={tilt_angle_mrad}"
    )
    echo(f"target defocus={target_defocus_um:.3f} um")
    echo(f"CTF X-tilt=({ctf_xtilt_x:.6f}, {ctf_xtilt_y:.6f})")
    echo(f"beam-tilt X-tilt=({beam_tilt_xtilt_x:.6f}, {beam_tilt_xtilt_y:.6f})")

    original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    echo(
        f"Saved XLensDeflector({xtilt_lens_index}): "
        f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
    )

    try:
        sem.GoToLowDoseArea("R")
        sem.SetImageShift(0, 0)

        echo("------------------------------------------------")
        echo("Beam-tilt autofocus")
        beam_defocus, raw = beam_tilt_autofocus(target_defocus_um)

        echo("------------------------------------------------")
        echo("CTF reference measurement")
        ctf_defocus, ctf_res = acquire_ctf_reference()

        physics = raw["legacy_defocus_um"]
        correction = raw["calibration_correction_um"]
        error_vs_target = beam_defocus - target_defocus_um
        delta_vs_ctf = beam_defocus - ctf_defocus
        physics_vs_ctf = physics - ctf_defocus

        echo("================================================")
        echo("FINAL RESULTS")
        echo(f"  target defocus:        {target_defocus_um:.4f} um")
        echo(f"  beam-tilt (calibrated): {beam_defocus:.4f} um")
        echo(f"  beam-tilt (physics):    {physics:.4f} um")
        echo(f"  calibration correction: {correction:.4f} um")
        echo(f"  Cs term:                {raw['cs_term_um']:.4f} um")
        echo(f"  CTF:                    {ctf_defocus:.4f} um ({ctf_res:.1f} A)")
        echo(f"  error vs target:        {error_vs_target:.4f} um")
        echo(f"  delta calibrated-CTF:   {delta_vs_ctf:.4f} um")
        echo(f"  delta physics-CTF:      {physics_vs_ctf:.4f} um")
        echo("================================================")

    finally:
        set_xtilt(float(original_xtilt[0]), float(original_xtilt[1]))
        echo(
            f"Restored XLensDeflector({xtilt_lens_index}) to "
            f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
        )

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
