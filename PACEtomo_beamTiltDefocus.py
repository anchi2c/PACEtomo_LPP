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

# Leave empty to use the legacy scalar beam-tilt equation. Set to a JSON file
# written by calibrate_beam_tilt_xtilt_matrix.py to use fitted coefficients.
calibration_file = ""

# Optional inline calibration. A calibration file, when set, overrides this.
calibration = {
    "model": "legacy_radial",
    "feature_names": [],
    "coefficients": [],
}

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
    """Return calibration from file, inline settings, or legacy fallback."""
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


def _legacy_defocus(raw, tilt_angle_mrad, legacy_divisor):
    shift_x = raw["shift_x_um"]
    shift_y = raw["shift_y_um"]
    displacement = np.sqrt(shift_x * shift_x + shift_y * shift_y)
    sign = 1.0 if shift_x == 0 else shift_x / abs(shift_x)
    return -1.0 * sign * displacement / (legacy_divisor * tilt_angle_mrad)


def defocus_from_raw(raw, tilt_angle_mrad=5.0, legacy_divisor=2.0,
                     calibration_data=None):
    """Convert raw beam-tilt shift diagnostics into defocus in microns."""
    calib = calibration_data if calibration_data is not None else get_calibration()
    if not calib or calib.get("model", "legacy_radial") == "legacy_radial":
        return _legacy_defocus(raw, tilt_angle_mrad, legacy_divisor)

    feature_names = calib.get("feature_names", [])
    coeffs = np.array(calib.get("coefficients", []), dtype=float)
    if len(feature_names) != len(coeffs):
        raise ValueError("Calibration feature_names and coefficients lengths differ")
    features = np.array([_feature_value(name, raw) for name in feature_names], dtype=float)
    return float(features.dot(coeffs))


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

    pixel_size_binned = float(_sem.ReportCurrentPixelSize("R"))
    binning = float(_sem.ReportBinning("R"))
    pixel_size_unbinned = pixel_size_binned / binning

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

    drift_x = disp_x2_px * pixel_size_unbinned
    drift_y = disp_y2_px * pixel_size_unbinned
    shift_x = (disp_x1_px - disp_x2_px / 2.0) * pixel_size_unbinned
    shift_y = (disp_y1_px - disp_y2_px / 2.0) * pixel_size_unbinned
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
        "pixel_size_unbinned_um": pixel_size_unbinned,
        "align_shift_1_x_px": disp_x1_px,
        "align_shift_1_y_px": disp_y1_px,
        "align_shift_2_x_px": disp_x2_px,
        "align_shift_2_y_px": disp_y2_px,
        "shift_x_um": shift_x,
        "shift_y_um": shift_y,
        "shift_abs_um": float(np.sqrt(shift_x * shift_x + shift_y * shift_y)),
        "drift_x_um": drift_x,
        "drift_y_um": drift_y,
        "drift_speed_x_nm_per_s": drift_x / elapsed if elapsed > 0 else 0.0,
        "drift_speed_y_nm_per_s": drift_y / elapsed if elapsed > 0 else 0.0,
        "elapsed_s": elapsed,
    }


def measure_defocus_with_diagnostics(tilt_angle_mrad=5.0,
                                     beam_tilt_correction=1.0,
                                     xtilt_x=None, xtilt_y=None,
                                     lens_index=2, beam_tilt_axis="x",
                                     legacy_divisor=2.0,
                                     calibration_data=None):
    """Measure defocus and return `(defocus, diagnostics)`."""
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
        defocus = defocus_from_raw(
            raw,
            tilt_angle_mrad=tilt_angle_mrad,
            legacy_divisor=legacy_divisor,
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
                    xtilt_x=None, xtilt_y=None, lens_index=2,
                    beam_tilt_axis="x", legacy_divisor=2.0,
                    calibration_data=None):
    """Measure defocus and return `(defocus, drift_speed_x, drift_speed_y)`."""
    defocus, raw = measure_defocus_with_diagnostics(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        xtilt_x=xtilt_x,
        xtilt_y=xtilt_y,
        lens_index=lens_index,
        beam_tilt_axis=beam_tilt_axis,
        legacy_divisor=legacy_divisor,
        calibration_data=calibration_data,
    )
    return defocus, raw["drift_speed_x_nm_per_s"], raw["drift_speed_y_nm_per_s"]


def autofocus_apply(target_defocus, cycles=2, tolerance_um=0.05,
                    tilt_angle_mrad=5.0, beam_tilt_correction=1.0,
                    xtilt_x=None, xtilt_y=None, lens_index=2,
                    beam_tilt_axis="x", legacy_divisor=2.0,
                    calibration_data=None):
    """Measure defocus by beam tilt and correct objective focus."""
    defocus = np.nan
    for cycle in range(1, int(cycles) + 1):
        defocus, speed_x, speed_y = measure_defocus(
            tilt_angle_mrad=tilt_angle_mrad,
            beam_tilt_correction=beam_tilt_correction,
            xtilt_x=xtilt_x,
            xtilt_y=xtilt_y,
            lens_index=lens_index,
            beam_tilt_axis=beam_tilt_axis,
            legacy_divisor=legacy_divisor,
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
