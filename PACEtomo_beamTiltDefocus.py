#!Python
# ===================================================================
# ScriptName     PACEtomo_beamTiltDefocus
# Purpose:       Shared calibrated beam-tilt defocus measurement helpers.
# ===================================================================

import json
import os

import numpy as np
import serialem as sem

############ SETTINGS ############

# Leave empty to use known physics (displacement/2beta - Cs*beta^2) only.
# Set to a JSON from calibrate_beam_tilt_xtilt_matrix.py to add a fitted
# residual correction on top of that physics base.
calibration_file = ""

# Optional inline calibration. A calibration file, when set, overrides this.
calibration = {
    "model": "physics_cs",
    "feature_names": [],
    "coefficients": [],
}

# Spherical aberration coefficient [mm]. Base beam-tilt equation:
#   defocus_um = -displacement_um / (2 * beta_rad) - Cs_um * beta_rad^2
#   beta_rad = tilt_angle_mrad * defocus_tilt_correction * 1e-3
# defocus_tilt_correction scales tilt in the defocus equation only (not SetBeamTilt).
spherical_aberration_mm = 2.7

########## END SETTINGS ##########

_sem = sem
_logger = None
_calibration_cache = None


def _echo(text):
    if _logger is not None:
        _logger(text)
    else:
        _sem.Echo(text)


def configure(sem_module=None, logger=None):
    """Set SerialEM module/logger from the importing script."""
    global _sem, _logger
    if sem_module is not None:
        _sem = sem_module
    _logger = logger


def load_calibration(path):
    """Load fitted calibration JSON."""
    with open(path, "r") as fh:
        return json.load(fh)


def get_calibration():
    """Return matrix calibration from file, inline settings, or physics-only default."""
    global _calibration_cache
    if calibration_file:
        if _calibration_cache is None:
            _calibration_cache = load_calibration(calibration_file)
        return _calibration_cache
    return calibration


def _feature_value(name, raw):
    """Feature values used by fitted calibration models."""
    shift_x = raw["shift_x_um"]
    shift_y = raw["shift_y_um"]
    xtilt_x = raw["xtilt_x"]
    xtilt_y = raw["xtilt_y"]
    beam_tilt_x0 = raw["beam_tilt_x0"]
    beam_tilt_y0 = raw["beam_tilt_y0"]
    tilt_step_x = raw["tilt_step_x"]
    tilt_step_y = raw["tilt_step_y"]
    values = {
        "1": 1.0,
        "shift_x_um": shift_x,
        "shift_y_um": shift_y,
        "shift_abs_um": raw["shift_abs_um"],
        "xtilt_x": xtilt_x,
        "xtilt_y": xtilt_y,
        "beam_tilt_x0": beam_tilt_x0,
        "beam_tilt_y0": beam_tilt_y0,
        "tilt_step_x": tilt_step_x,
        "tilt_step_y": tilt_step_y,
        "shift_x_um*xtilt_x": shift_x * xtilt_x,
        "shift_x_um*xtilt_y": shift_x * xtilt_y,
        "shift_y_um*xtilt_x": shift_y * xtilt_x,
        "shift_y_um*xtilt_y": shift_y * xtilt_y,
        "shift_x_um*beam_tilt_x0": shift_x * beam_tilt_x0,
        "shift_x_um*beam_tilt_y0": shift_x * beam_tilt_y0,
        "shift_y_um*beam_tilt_x0": shift_y * beam_tilt_x0,
        "shift_y_um*beam_tilt_y0": shift_y * beam_tilt_y0,
    }
    if name not in values:
        raise KeyError(f"Unknown calibration feature '{name}'")
    return float(values[name])


def _signed_displacement_um(raw, beam_tilt_axis="x"):
    """Signed image displacement along the beam-tilt measurement axis [um]."""
    if beam_tilt_axis.lower() == "y":
        return float(raw["shift_y_um"])
    return float(raw["shift_x_um"])


def _physics_cs_defocus(raw, tilt_angle_mrad, cs_mm, beam_tilt_axis="x",
                        defocus_tilt_correction=1.0):
    """
    Known physics beam-tilt defocus with Cs [um].

    defocus_um = -displacement_um / (2 * beta_rad) - Cs_um * beta_rad^2
    beta_rad = tilt_angle_mrad * defocus_tilt_correction * 1e-3

    defocus_tilt_correction enters the defocus equation only. SetBeamTilt uses
    beam_tilt_correction separately (often ~1.73); many scopes need ~0.95 here.
    """
    beta_rad = float(tilt_angle_mrad) * float(defocus_tilt_correction) * 1e-3
    signed_displacement = _signed_displacement_um(raw, beam_tilt_axis)
    cs_um = float(cs_mm) * 1000.0
    linear_term_um = -signed_displacement / (2.0 * beta_rad)
    cs_term_um = cs_um * beta_rad * beta_rad
    defocus_um = linear_term_um - cs_term_um
    return defocus_um, linear_term_um, cs_term_um, beta_rad


def _fitted_correction(raw, calib):
    """Empirical residual from matrix calibration (added on top of physics base)."""
    feature_names = calib.get("feature_names", [])
    coeffs = np.array(calib.get("coefficients", []), dtype=float)
    if len(feature_names) != len(coeffs):
        raise ValueError("Calibration feature_names and coefficients lengths differ")
    features = np.array([_feature_value(name, raw) for name in feature_names], dtype=float)
    return float(features.dot(coeffs))


def defocus_from_raw(raw, tilt_angle_mrad=5.0, cs_mm=None, beam_tilt_axis="x",
                     defocus_tilt_correction=None, calibration_data=None):
    """
    Convert raw beam-tilt diagnostics to defocus [um].

    Always starts from the physics+Cs base. A loaded matrix calibration adds a
    fitted residual on top; it does not replace the known equation.
    """
    if cs_mm is None:
        cs_mm = spherical_aberration_mm
    calib = calibration_data if calibration_data is not None else get_calibration()
    if defocus_tilt_correction is None:
        if calib and "defocus_tilt_correction" in calib:
            defocus_tilt_correction = float(calib["defocus_tilt_correction"])
        else:
            defocus_tilt_correction = float(raw.get("defocus_tilt_correction", 1.0))
    base_um, _, _, _ = _physics_cs_defocus(
        raw, tilt_angle_mrad, cs_mm, beam_tilt_axis=beam_tilt_axis,
        defocus_tilt_correction=defocus_tilt_correction,
    )
    if not calib:
        return float(base_um)
    model = calib.get("model", "physics_cs")
    if model == "physics_cs":
        return float(base_um)
    if model == "linear_xtilt_residual":
        return float(base_um + _fitted_correction(raw, calib))
    if model == "linear_xtilt_beam_tilt":
        # Older absolute fit (replaces base); kept for existing JSON files.
        return float(_fitted_correction(raw, calib))
    return float(base_um)


def legacy_physics_diagnostics(raw, tilt_angle_mrad=5.0, cs_mm=None,
                               beam_tilt_axis="x",
                               defocus_tilt_correction=1.0):
    """Term breakdown for the physics+Cs equation."""
    if cs_mm is None:
        cs_mm = spherical_aberration_mm
    defocus_um, linear_term_um, cs_term_um, beta_rad = _physics_cs_defocus(
        raw, tilt_angle_mrad, cs_mm, beam_tilt_axis=beam_tilt_axis,
        defocus_tilt_correction=defocus_tilt_correction,
    )
    return {
        "beta_rad": float(beta_rad),
        "defocus_tilt_correction": float(defocus_tilt_correction),
        "linear_term_um": float(linear_term_um),
        "cs_term_um": float(cs_term_um),
        "legacy_defocus_um": float(defocus_um),
        "spherical_aberration_mm": float(cs_mm),
    }


def empirical_physics_cs_defocus(raw, tilt_angle_mrad=5.0, cs_mm=None,
                                 beam_tilt_axis="x",
                                 defocus_tilt_correction=1.0):
    """
    Same physics+Cs equation with a trial defocus_tilt_correction in beta.

    Typically ~0.95 instead of beam_tilt_correction (~1.73) on SetBeamTilt.
    """
    if cs_mm is None:
        cs_mm = spherical_aberration_mm
    defocus_um, _, _, _ = _physics_cs_defocus(
        raw, tilt_angle_mrad, cs_mm, beam_tilt_axis=beam_tilt_axis,
        defocus_tilt_correction=defocus_tilt_correction,
    )
    return float(defocus_um)


# physics_cs_diagnostics kept as alias for older call sites
physics_cs_diagnostics = legacy_physics_diagnostics


def calibration_range_warnings(raw, calibration_data=None):
    """Return warnings when raw state is outside the fitted calibration range."""
    calib = calibration_data if calibration_data is not None else get_calibration()
    if not calib:
        return []
    warnings = []
    for name, limits in calib.get("raw_ranges", {}).items():
        if name not in raw or len(limits) != 2:
            continue
        value = float(raw[name])
        lo = float(limits[0])
        hi = float(limits[1])
        if value < lo or value > hi:
            warnings.append(
                f"{name}={value:.6g} outside calibration range [{lo:.6g}, {hi:.6g}]"
            )
    return warnings


def measure_raw(tilt_angle_mrad=5.0, beam_tilt_correction=1.0,
                beam_tilt_axis="x"):
    """
    Acquire +tilt, -tilt, +tilt images and return raw signed shift diagnostics.

    The third image estimates drift during the pair, matching the existing
    PACEtomo beam-tilt autofocus sequence.
    """
    beam_tilt = _sem.ReportBeamTilt()
    tilt_x_orig = float(beam_tilt[0])
    tilt_y_orig = float(beam_tilt[1])
    tilt_step = float(beam_tilt_correction) * float(tilt_angle_mrad)

    if beam_tilt_axis.lower() == "y":
        tilt_x_plus = tilt_x_orig
        tilt_x_minus = tilt_x_orig
        tilt_y_plus = tilt_y_orig + tilt_step
        tilt_y_minus = tilt_y_orig - tilt_step
        tilt_step_x = 0.0
        tilt_step_y = tilt_step
    else:
        tilt_x_plus = tilt_x_orig + tilt_step
        tilt_x_minus = tilt_x_orig - tilt_step
        tilt_y_plus = tilt_y_orig
        tilt_y_minus = tilt_y_orig
        tilt_step_x = tilt_step
        tilt_step_y = 0.0

    focus_camera = "F"
    pixel_size_binned_nm = float(_sem.ReportCurrentPixelSize(focus_camera))
    focus_binning = float(_sem.ReportBinning(focus_camera))
    # SerialEM returns pixel size in nm; defocus equation uses displacement [um].
    pixel_size_unbinned_nm = pixel_size_binned_nm / focus_binning

    _sem.SetBeamTilt(tilt_x_plus, tilt_y_plus)
    _sem.F()
    _sem.ResetClock()
    _sem.Copy("A", "L")

    _sem.SetBeamTilt(tilt_x_minus, tilt_y_minus)
    _sem.F()
    _sem.AlignTo("L", 1)
    align_shift_1 = _sem.ReportAlignShift()
    disp_x1_px = float(align_shift_1[0])
    disp_y1_px = float(align_shift_1[1])

    _sem.SetBeamTilt(tilt_x_plus, tilt_y_plus)
    _sem.F()
    elapsed = float(_sem.ReportClock())

    _sem.SetBeamTilt(tilt_x_orig, tilt_y_orig)
    _sem.AlignTo("L", 1)
    align_shift_2 = _sem.ReportAlignShift()
    disp_x2_px = float(align_shift_2[0])
    disp_y2_px = float(align_shift_2[1])

    drift_x_nm = disp_x2_px * pixel_size_unbinned_nm
    drift_y_nm = disp_y2_px * pixel_size_unbinned_nm
    shift_x_um = (disp_x1_px - disp_x2_px / 2.0) * pixel_size_unbinned_nm / 1000.0
    shift_y_um = (disp_y1_px - disp_y2_px / 2.0) * pixel_size_unbinned_nm / 1000.0
    pixel_size_um = pixel_size_unbinned_nm / 1000.0
    if beam_tilt_axis.lower() == "y":
        minus_vs_plus_ref_um = disp_y1_px * pixel_size_um
        drift_return_um = disp_y2_px * pixel_size_um
        shift_axis_um = shift_y_um
    else:
        minus_vs_plus_ref_um = disp_x1_px * pixel_size_um
        drift_return_um = disp_x2_px * pixel_size_um
        shift_axis_um = shift_x_um
    # Drift-corrected +/- branch shifts about zero beam tilt (symmetric pair).
    plus_branch_shift_um = 0.5 * shift_axis_um
    minus_branch_shift_um = -0.5 * shift_axis_um
    xtilt = _sem.ReportXLensDeflector(2)

    return {
        "beam_tilt_x0": tilt_x_orig,
        "beam_tilt_y0": tilt_y_orig,
        "tilt_step_x": tilt_step_x,
        "tilt_step_y": tilt_step_y,
        "tilt_angle_mrad": float(tilt_angle_mrad),
        "beam_tilt_correction": float(beam_tilt_correction),
        "xtilt_x": float(xtilt[0]),
        "xtilt_y": float(xtilt[1]),
        "focus_camera": focus_camera,
        "focus_binning": focus_binning,
        "focus_pixel_size_binned_nm": pixel_size_binned_nm,
        "focus_pixel_size_unbinned_nm": pixel_size_unbinned_nm,
        "align_shift_1_x_px": disp_x1_px,
        "align_shift_1_y_px": disp_y1_px,
        "align_shift_2_x_px": disp_x2_px,
        "align_shift_2_y_px": disp_y2_px,
        "shift_x_um": shift_x_um,
        "shift_y_um": shift_y_um,
        "shift_axis_um": float(shift_axis_um),
        "plus_branch_shift_um": float(plus_branch_shift_um),
        "minus_branch_shift_um": float(minus_branch_shift_um),
        "minus_vs_plus_ref_um": float(minus_vs_plus_ref_um),
        "drift_return_um": float(drift_return_um),
        "shift_abs_um": float(np.sqrt(shift_x_um * shift_x_um + shift_y_um * shift_y_um)),
        "drift_x_nm": drift_x_nm,
        "drift_y_nm": drift_y_nm,
        "drift_speed_x_nm_per_s": drift_x_nm / elapsed if elapsed > 0 else 0.0,
        "drift_speed_y_nm_per_s": drift_y_nm / elapsed if elapsed > 0 else 0.0,
        "elapsed_s": elapsed,
    }


def measure_defocus_with_diagnostics(tilt_angle_mrad=5.0,
                                     beam_tilt_correction=1.0,
                                     defocus_tilt_correction=None,
                                     xtilt_x=None, xtilt_y=None,
                                     lens_index=2, beam_tilt_axis="x",
                                     cs_mm=None, calibration_data=None):
    """Measure defocus and return `(defocus, diagnostics)`."""
    if cs_mm is None:
        cs_mm = spherical_aberration_mm
    original_xtilt = _sem.ReportXLensDeflector(lens_index)
    original_beam_tilt = _sem.ReportBeamTilt()
    try:
        if xtilt_x is not None and xtilt_y is not None:
            _sem.SetXLensDeflector(lens_index, float(xtilt_x), float(xtilt_y))
        raw = measure_raw(
            tilt_angle_mrad=tilt_angle_mrad,
            beam_tilt_correction=beam_tilt_correction,
            beam_tilt_axis=beam_tilt_axis,
        )
        if defocus_tilt_correction is not None:
            raw["defocus_tilt_correction"] = float(defocus_tilt_correction)
        raw.update(legacy_physics_diagnostics(
            raw,
            tilt_angle_mrad=tilt_angle_mrad,
            cs_mm=cs_mm,
            beam_tilt_axis=beam_tilt_axis,
            defocus_tilt_correction=(
                defocus_tilt_correction
                if defocus_tilt_correction is not None
                else float(raw.get("defocus_tilt_correction", 1.0))
            ),
        ))
        calib = calibration_data if calibration_data is not None else get_calibration()
        model = calib.get("model", "physics_cs") if calib else "physics_cs"
        correction_um = 0.0
        if model == "linear_xtilt_residual":
            correction_um = _fitted_correction(raw, calib)
        raw["calibration_correction_um"] = float(correction_um)
        defocus = defocus_from_raw(
            raw,
            tilt_angle_mrad=tilt_angle_mrad,
            cs_mm=cs_mm,
            beam_tilt_axis=beam_tilt_axis,
            defocus_tilt_correction=defocus_tilt_correction,
            calibration_data=calibration_data,
        )
        raw["defocus_um"] = float(defocus)
        raw["calibration_warnings"] = calibration_range_warnings(
            raw, calibration_data=calibration_data
        )
        for warning in raw["calibration_warnings"]:
            _echo(f"WARNING: Beam-tilt calibration: {warning}")
        return float(defocus), raw
    finally:
        _sem.SetBeamTilt(float(original_beam_tilt[0]), float(original_beam_tilt[1]))
        _sem.SetXLensDeflector(lens_index, float(original_xtilt[0]), float(original_xtilt[1]))


def measure_defocus(tilt_angle_mrad=5.0, beam_tilt_correction=1.0,
                    defocus_tilt_correction=None,
                    xtilt_x=None, xtilt_y=None, lens_index=2,
                    beam_tilt_axis="x", cs_mm=None,
                    calibration_data=None):
    """Measure defocus and return `(defocus, drift_speed_x, drift_speed_y)`."""
    defocus, raw = measure_defocus_with_diagnostics(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        defocus_tilt_correction=defocus_tilt_correction,
        xtilt_x=xtilt_x,
        xtilt_y=xtilt_y,
        lens_index=lens_index,
        beam_tilt_axis=beam_tilt_axis,
        cs_mm=cs_mm,
        calibration_data=calibration_data,
    )
    return defocus, raw["drift_speed_x_nm_per_s"], raw["drift_speed_y_nm_per_s"]


def autofocus_apply(target_defocus, cycles=2, tolerance_um=0.05,
                    tilt_angle_mrad=5.0, beam_tilt_correction=1.0,
                    defocus_tilt_correction=None,
                    xtilt_x=None, xtilt_y=None, lens_index=2,
                    beam_tilt_axis="x", cs_mm=None,
                    calibration_data=None):
    """Measure defocus by beam tilt and correct objective focus."""
    defocus = np.nan
    for cycle in range(1, int(cycles) + 1):
        defocus, speed_x, speed_y = measure_defocus(
            tilt_angle_mrad=tilt_angle_mrad,
            beam_tilt_correction=beam_tilt_correction,
            defocus_tilt_correction=defocus_tilt_correction,
            xtilt_x=xtilt_x,
            xtilt_y=xtilt_y,
            lens_index=lens_index,
            beam_tilt_axis=beam_tilt_axis,
            cs_mm=cs_mm,
            calibration_data=calibration_data,
        )
        error = float(target_defocus) - defocus
        _echo(
            f"Autofocus {cycle}/{int(cycles)}: measured={defocus:.4f} um, "
            f"target={float(target_defocus):.3f} um, error={error:.3f} um, "
            f"drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
        )
        if abs(error) <= float(tolerance_um):
            return defocus
        _sem.ChangeFocus(error)
    return defocus


def calibration_path_in_current_dir(filename):
    """Convenience for SerialEM working directories."""
    if not filename:
        return ""
    if os.path.isabs(filename):
        return filename
    return os.path.abspath(filename)
