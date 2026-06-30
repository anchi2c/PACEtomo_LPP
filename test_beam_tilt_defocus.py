#!Python
# ===================================================================
# ScriptName     test_beam_tilt_defocus
# Purpose:       Validate PACEtomo_beamTiltDefocus vs CTF using the same
#                parameters as calibrate_beam_tilt_scaling.py.
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

target_defocus_values = [-1.0, -2.0, -3.0, -4.0, -5.0]

ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2
target_defocus_tolerance_um = 0.05
max_defocus_adjust_iterations = 5

xtilt_lens_index = 2
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867
beam_tilt_xtilt_x = 0.0
beam_tilt_xtilt_y = 0.0

# Point at scaling JSON to test calibration; leave "" for physics-only.
calibration_file = r""

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


def set_target_defocus(target_defocus):
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    current_defocus = float("nan")
    for attempt in range(1, max_defocus_adjust_iterations + 1):
        current_defocus, _ = acquire_ctf_reference()
        error = target_defocus - current_defocus
        echo(
            f"Target defocus {attempt}/{max_defocus_adjust_iterations}: "
            f"current={current_defocus:.3f} um, target={target_defocus:.3f} um"
        )
        if abs(error) <= target_defocus_tolerance_um:
            return current_defocus
        sem.ChangeFocus(error)
    echo("WARNING: Target defocus not reached within tolerance.")
    return current_defocus


def main():
    sem.SuppressReports()
    echo("##### test_beam_tilt_defocus #####")
    echo(f"calibration_file={btdef.calibration_file or '(none)'}")
    echo(
        f"beam_tilt_correction={beam_tilt_correction}, "
        f"defocus_tilt_correction={defocus_tilt_correction}, "
        f"tilt_angle_mrad={tilt_angle_mrad}"
    )
    echo(f"CTF X-tilt=({ctf_xtilt_x:.6f}, {ctf_xtilt_y:.6f})")
    echo(f"beam-tilt X-tilt=({beam_tilt_xtilt_x:.6f}, {beam_tilt_xtilt_y:.6f})")

    original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    echo(
        f"Saved XLensDeflector({xtilt_lens_index}): "
        f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
    )

    deltas_physics = []
    deltas_calibrated = []
    try:
        for target_defocus in target_defocus_values:
            echo("------------------------------------------------")
            echo(f"Target defocus {float(target_defocus):.3f} um")
            reached = set_target_defocus(target_defocus)
            echo(f"CTF after focus adjust: {reached:.3f} um")

            defocus, raw = btdef.measure_defocus_with_diagnostics(
                tilt_angle_mrad=tilt_angle_mrad,
                beam_tilt_correction=beam_tilt_correction,
                defocus_tilt_correction=defocus_tilt_correction,
                xtilt_x=beam_tilt_xtilt_x,
                xtilt_y=beam_tilt_xtilt_y,
                lens_index=xtilt_lens_index,
                beam_tilt_axis=beam_tilt_axis,
            )
            ctf_defocus, ctf_res = acquire_ctf_reference()

            physics = raw["legacy_defocus_um"]
            correction = raw["calibration_correction_um"]
            delta_physics = physics - ctf_defocus
            delta_calibrated = defocus - ctf_defocus
            deltas_physics.append(delta_physics)
            deltas_calibrated.append(delta_calibrated)

            echo(
                f"physics={physics:.4f}, Cs={raw['cs_term_um']:.4f}, "
                f"correction={correction:.4f}, calibrated={defocus:.4f}, "
                f"CTF={ctf_defocus:.4f} ({ctf_res:.1f} A), "
                f"delta(physics)={delta_physics:.4f}, "
                f"delta(calibrated)={delta_calibrated:.4f}, "
                f"drift=({raw['drift_speed_x_nm_per_s']:.2f}, "
                f"{raw['drift_speed_y_nm_per_s']:.2f}) nm/s"
            )

        echo("================================================")
        echo(
            f"mean delta(physics)={sum(deltas_physics) / len(deltas_physics):.4f} um, "
            f"mean delta(calibrated)={sum(deltas_calibrated) / len(deltas_calibrated):.4f} um"
        )
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
