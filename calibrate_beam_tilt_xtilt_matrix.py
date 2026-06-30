#!Python
# ===================================================================
# ScriptName     calibrate_beam_tilt_xtilt_matrix
# Purpose:       Measure how beam tilt, residual beam tilt, and
#                XLensDeflector(2) affect beam-tilt defocus autofocus.
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import csv
import json
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

import PACEtomo_beamTiltDefocus as btdef

############ SETTINGS ############

# Beam tilt used for the autofocus-like measurement [mrad].
tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
beam_tilt_axis = "x"

# Spherical aberration [mm] for defocus = displacement/(2*beta) - Cs*beta^2
spherical_aberration_mm = 2.7

# Empirical beta scale in physics+Cs:
#   defocus = displacement/(2*beta*f) - Cs*(beta*f)^2  (beta = tilt_angle_mrad * 1e-3)
# Separate from hardware beam_tilt_correction (SetBeamTilt scale).
empirical_defocus_factor = 1.0

# Optional sweep grid for beta-factor fit. Leave empty for 0.7-1.3.
empirical_defocus_factors_for_fit = []

# Nominal X-tilt used as the center of the beam-tilt autofocus sweep.
nominal_xtilt_x = 0.0
nominal_xtilt_y = 0.0

# X-tilt used for CTF reference defocus. This can be the high/known-good
# PACEtomo.py ctfXtiltX/Y setting, independent of the X-tilt being tested.
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867

# Test X-tilt around nominal. Start one axis at a time; enable the 2D grid
# after the one-axis tests show the sweep range is safe.
xtilt_step = 0.0005
xtilt_offsets_x = [-2, -1, 0, 1, 2]
xtilt_offsets_y = [0]
include_xtilt_2d_grid = False

# Residual beam tilt/coma offsets around the current microscope setting.
residual_beam_tilt_step = 0.25
#beam_tilt_offsets_x = [-1, 0, 1]
beam_tilt_offsets_x = [0]
beam_tilt_offsets_y = [0]
include_beam_tilt_2d_grid = False

# Reference defocus values [microns].
target_defocus_values = [-1.0, -2.0, -3.0]

# Optional stage Z offsets [microns]. Leave at [0.0] unless you want to test
# height-induced defocus separately from objective-lens defocus.
stage_z_offsets_um = [0.0]
use_stage_z_offsets = False
stage_z_backlash_um = 1.0

# CTF search range [microns].
ctf_defocus_lo = -12.0
ctf_defocus_hi = -0.2
ctf_resolution_max_A = 20.0

# Defocus convergence settings for setting each target defocus.
target_defocus_tolerance_um = 0.025
max_defocus_adjust_iterations = 5

# Repeats per condition.
repeats = 1

# Reject rows from fitting if drift exceeds this speed. Set <= 0 to disable.
max_drift_nm_per_s = 2.0

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r"X:\k3f_leginonframes\p26jun29a\xtilt_calib_test"

# Output file names (written under save_dir when set).
csv_measurements = "beam_tilt_xtilt_matrix_measurements.csv"
csv_summary = "beam_tilt_xtilt_matrix_summary.csv"
calibration_json = "beam_tilt_xtilt_matrix_calibration.json"
plot_file = "beam_tilt_xtilt_matrix_plots.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


btdef.configure(sem_module=sem, logger=echo)


def prepare_output_paths():
    """Resolve output paths under save_dir and create the directory if needed."""
    global csv_measurements, csv_summary, calibration_json, plot_file
    if not save_dir:
        return
    os.makedirs(save_dir, exist_ok=True)
    csv_measurements = os.path.join(save_dir, csv_measurements)
    csv_summary = os.path.join(save_dir, csv_summary)
    calibration_json = os.path.join(save_dir, calibration_json)
    plot_file = os.path.join(save_dir, plot_file)


def legacy_ctf_ratio(legacy_defocus, ctf_defocus):
    """Ratio legacy / CTF defocus; NaN when CTF is zero or values are invalid."""
    if (not np.isfinite(legacy_defocus) or not np.isfinite(ctf_defocus)
            or abs(float(ctf_defocus)) < 1e-6):
        return np.nan
    return float(legacy_defocus) / float(ctf_defocus)


def row_to_raw(row):
    return {
        "shift_x_um": float(row["shift_x_um"]),
        "shift_y_um": float(row["shift_y_um"]),
        "shift_abs_um": float(row["shift_abs_um"]),
        "xtilt_x": float(row["xtilt_x"]),
        "xtilt_y": float(row["xtilt_y"]),
        "beam_tilt_x0": float(row["beam_tilt_x0"]),
        "beam_tilt_y0": float(row["beam_tilt_y0"]),
        "tilt_step_x": float(row["tilt_step_x"]),
        "tilt_step_y": float(row["tilt_step_y"]),
    }


def fit_empirical_beta_factor_zero_crossing(fit_rows, factors):
    """Fit beta factor where mean(empirical physics+Cs - CTF) crosses zero."""
    if len(fit_rows) == 0 or len(factors) < 2:
        return np.nan, np.nan, np.nan
    mean_deltas = []
    for factor in factors:
        deltas = []
        for row in fit_rows:
            empirical = btdef.empirical_physics_cs_defocus(
                row_to_raw(row),
                tilt_angle_mrad=tilt_angle_mrad,
                cs_mm=spherical_aberration_mm,
                beam_tilt_axis=beam_tilt_axis,
                beta_factor=factor,
            )
            deltas.append(empirical - float(row["ctf_defocus_um"]))
        mean_deltas.append(float(np.mean(deltas)))
    slope, intercept = np.polyfit(np.array(factors, dtype=float),
                                  np.array(mean_deltas, dtype=float), 1)
    if slope == 0:
        return slope, intercept, np.nan
    return slope, intercept, -intercept / slope


def fit_empirical_beta_factor_min_rms(fit_rows, factors):
    """Grid search for beta factor minimizing RMS vs CTF (nonlinear with Cs)."""
    if len(fit_rows) == 0 or len(factors) == 0:
        return np.nan, np.nan
    best_factor = np.nan
    best_rms = np.nan
    for factor in factors:
        rms = empirical_rms_error_um(fit_rows, factor)
        if np.isfinite(rms) and (not np.isfinite(best_rms) or rms < best_rms):
            best_factor = float(factor)
            best_rms = float(rms)
    return best_factor, best_rms


def apply_empirical_fitted_columns(rows, fitted_factor):
    """Add recalculated physics+Cs defocus using the fitted beta factor."""
    for row in rows:
        empirical_fitted = btdef.empirical_physics_cs_defocus(
            row_to_raw(row),
            tilt_angle_mrad=tilt_angle_mrad,
            cs_mm=spherical_aberration_mm,
            beam_tilt_axis=beam_tilt_axis,
            beta_factor=fitted_factor,
        )
        ctf = float(row["ctf_defocus_um"])
        row["fitted_empirical_defocus_factor"] = float(fitted_factor)
        row["empirical_fitted_defocus_um"] = float(empirical_fitted)
        row["delta_empirical_fitted_minus_ctf_um"] = empirical_fitted - ctf
        row["ratio_empirical_fitted_over_ctf"] = legacy_ctf_ratio(
            empirical_fitted, ctf
        )


def axis_or_grid(x_steps, y_steps, step_size, include_grid):
    if include_grid:
        return [(x * step_size, y * step_size) for x in x_steps for y in y_steps]
    offsets = [(x * step_size, 0.0) for x in x_steps]
    offsets.extend((0.0, y * step_size) for y in y_steps if y != 0)
    seen = set()
    unique = []
    for x, y in offsets:
        key = (round(x, 12), round(y, 12))
        if key not in seen:
            seen.add(key)
            unique.append((x, y))
    return unique


def run_ctffind():
    sem.NoMessageBoxOnError(1)
    try:
        cfind = sem.CtfFind("A", ctf_defocus_lo, ctf_defocus_hi)
    finally:
        sem.NoMessageBoxOnError(0)
    if len(cfind) == 0:
        return np.nan, np.nan
    return float(cfind[0]), float(cfind[-1])


def acquire_ctf_reference():
    """Acquire CTF reference defocus at the dedicated high-X-tilt setting."""
    xtilt = sem.ReportXLensDeflector(2)
    try:
        sem.SetXLensDeflector(2, ctf_xtilt_x, ctf_xtilt_y)
        sem.L()
        return run_ctffind()
    finally:
        sem.SetXLensDeflector(2, float(xtilt[0]), float(xtilt[1]))


def set_target_defocus(target_defocus):
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    current_defocus = np.nan
    for attempt in range(1, max_defocus_adjust_iterations + 1):
        current_defocus, ctf_res = acquire_ctf_reference()
        if not np.isfinite(current_defocus):
            sem.Echo("WARNING: CtfFind failed while setting target defocus.")
            continue
        error = float(target_defocus) - current_defocus
        sem.Echo(
            f"Target defocus check {attempt}/{max_defocus_adjust_iterations}: "
            f"current={current_defocus:.3f} um, target={float(target_defocus):.3f} um, "
            f"error={error:.3f} um, res={ctf_res:.2f} A"
        )
        if abs(error) <= target_defocus_tolerance_um:
            return current_defocus
        sem.ChangeFocus(error)
    sem.Echo("WARNING: Target defocus not reached within tolerance.")
    return current_defocus


def move_stage_z(original_stage, z_offset_um):
    if not use_stage_z_offsets or abs(float(z_offset_um)) == 0:
        return
    # SerialEM stage Z is reported in microns on current supported versions.
    target_z = float(original_stage[2]) + float(z_offset_um)
    backlash = abs(float(stage_z_backlash_um))
    if backlash > 0:
        direction = 1.0 if float(z_offset_um) >= 0 else -1.0
        sem.MoveStageTo(float(original_stage[0]), float(original_stage[1]),
                        target_z + direction * backlash)
    sem.MoveStageTo(float(original_stage[0]), float(original_stage[1]), target_z)


def restore_stage(original_stage):
    if not use_stage_z_offsets:
        return
    target_z = float(original_stage[2])
    backlash = abs(float(stage_z_backlash_um))
    if backlash > 0:
        sem.MoveStageTo(float(original_stage[0]), float(original_stage[1]),
                        target_z - backlash)
    sem.MoveStageTo(float(original_stage[0]), float(original_stage[1]), target_z)


def measurement_fields():
    return [
        "timestamp",
        "target_defocus_um",
        "stage_z_offset_um",
        "repeat",
        "xtilt_offset_x",
        "xtilt_offset_y",
        "xtilt_x",
        "xtilt_y",
        "ctf_xtilt_x",
        "ctf_xtilt_y",
        "beam_tilt_offset_x",
        "beam_tilt_offset_y",
        "beam_tilt_x0",
        "beam_tilt_y0",
        "tilt_angle_mrad",
        "beam_tilt_correction",
        "tilt_step_x",
        "tilt_step_y",
        "focus_binning",
        "focus_pixel_size_binned_nm",
        "focus_pixel_size_unbinned_nm",
        "shift_x_um",
        "shift_y_um",
        "shift_abs_um",
        "beta_rad",
        "linear_term_um",
        "cs_term_um",
        "spherical_aberration_mm",
        "legacy_defocus_um",
        "empirical_defocus_factor",
        "empirical_defocus_um",
        "delta_empirical_minus_ctf_um",
        "ratio_empirical_over_ctf",
        "fitted_empirical_defocus_factor",
        "empirical_fitted_defocus_um",
        "delta_empirical_fitted_minus_ctf_um",
        "ratio_empirical_fitted_over_ctf",
        "ctf_defocus_um",
        "delta_legacy_minus_ctf_um",
        "ratio_legacy_over_ctf",
        "ctf_resolution_A",
        "drift_speed_x_nm_per_s",
        "drift_speed_y_nm_per_s",
        "drift_speed_abs_nm_per_s",
        "elapsed_s",
        "fit_used",
    ]


def append_measurement(path, row):
    exists = False
    try:
        with open(path, "r"):
            exists = True
    except FileNotFoundError:
        pass
    with open(path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=measurement_fields())
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def initialize_measurements(path):
    """Start a fresh raw CSV so field changes cannot corrupt the header."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=measurement_fields())
        writer.writeheader()


def write_measurements_csv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=measurement_fields())
        writer.writeheader()
        writer.writerows(rows)


def collect_measurements():
    rows = []
    original_xtilt = sem.ReportXLensDeflector(2)
    original_beam_tilt = sem.ReportBeamTilt()
    original_stage = sem.ReportStageXYZ()
    xtilt_offsets = axis_or_grid(
        xtilt_offsets_x, xtilt_offsets_y, xtilt_step, include_xtilt_2d_grid
    )
    beam_offsets = axis_or_grid(
        beam_tilt_offsets_x,
        beam_tilt_offsets_y,
        residual_beam_tilt_step,
        include_beam_tilt_2d_grid,
    )

    sem.Echo(f"X-tilt offsets: {xtilt_offsets}")
    sem.Echo(f"Residual beam-tilt offsets: {beam_offsets}")

    try:
        for stage_z_offset in stage_z_offsets_um:
            restore_stage(original_stage)
            move_stage_z(original_stage, stage_z_offset)
            for target_defocus in target_defocus_values:
                sem.Echo("------------------------------------------------")
                sem.Echo(
                    f"Target defocus {float(target_defocus):.3f} um, "
                    f"stage_z_offset={float(stage_z_offset):.3f} um"
                )
                reached = set_target_defocus(target_defocus)
                sem.Echo(f"Defocus after adjustment loop: {reached:.3f} um")

                for xt_off_x, xt_off_y in xtilt_offsets:
                    xtilt_x = nominal_xtilt_x + xt_off_x
                    xtilt_y = nominal_xtilt_y + xt_off_y
                    sem.SetXLensDeflector(2, xtilt_x, xtilt_y)

                    for bt_off_x, bt_off_y in beam_offsets:
                        for repeat in range(1, int(repeats) + 1):
                            sem.SetBeamTilt(
                                float(original_beam_tilt[0]) + bt_off_x,
                                float(original_beam_tilt[1]) + bt_off_y,
                            )
                            raw = btdef.measure_raw(
                                tilt_angle_mrad=tilt_angle_mrad,
                                beam_tilt_correction=beam_tilt_correction,
                                beam_tilt_axis=beam_tilt_axis,
                            )
                            legacy = btdef.legacy_physics_diagnostics(
                                raw,
                                tilt_angle_mrad=tilt_angle_mrad,
                                cs_mm=spherical_aberration_mm,
                                beam_tilt_axis=beam_tilt_axis,
                            )
                            legacy_defocus = legacy["legacy_defocus_um"]
                            empirical_defocus = btdef.empirical_physics_cs_defocus(
                                raw,
                                tilt_angle_mrad=tilt_angle_mrad,
                                cs_mm=spherical_aberration_mm,
                                beam_tilt_axis=beam_tilt_axis,
                                beta_factor=empirical_defocus_factor,
                            )
                            ctf_defocus, ctf_res = acquire_ctf_reference()
                            delta_legacy_ctf = legacy_defocus - ctf_defocus
                            ratio_legacy_ctf = legacy_ctf_ratio(
                                legacy_defocus, ctf_defocus
                            )
                            delta_empirical_ctf = empirical_defocus - ctf_defocus
                            ratio_empirical_ctf = legacy_ctf_ratio(
                                empirical_defocus, ctf_defocus
                            )
                            drift_abs = float(np.hypot(
                                raw["drift_speed_x_nm_per_s"],
                                raw["drift_speed_y_nm_per_s"],
                            ))
                            fit_used = (
                                np.isfinite(ctf_defocus)
                                and np.isfinite(ctf_res)
                                and ctf_res <= ctf_resolution_max_A
                                and (max_drift_nm_per_s <= 0
                                     or drift_abs <= max_drift_nm_per_s)
                            )
                            row = {
                                "timestamp": datetime.now().isoformat(timespec="seconds"),
                                "target_defocus_um": float(target_defocus),
                                "stage_z_offset_um": float(stage_z_offset),
                                "repeat": int(repeat),
                                "xtilt_offset_x": xt_off_x,
                                "xtilt_offset_y": xt_off_y,
                                "xtilt_x": xtilt_x,
                                "xtilt_y": xtilt_y,
                                "ctf_xtilt_x": ctf_xtilt_x,
                                "ctf_xtilt_y": ctf_xtilt_y,
                                "beam_tilt_offset_x": bt_off_x,
                                "beam_tilt_offset_y": bt_off_y,
                                "beam_tilt_x0": raw["beam_tilt_x0"],
                                "beam_tilt_y0": raw["beam_tilt_y0"],
                                "tilt_angle_mrad": raw["tilt_angle_mrad"],
                                "beam_tilt_correction": raw["beam_tilt_correction"],
                                "tilt_step_x": raw["tilt_step_x"],
                                "tilt_step_y": raw["tilt_step_y"],
                                "focus_binning": raw["focus_binning"],
                                "focus_pixel_size_binned_nm": raw["focus_pixel_size_binned_nm"],
                                "focus_pixel_size_unbinned_nm": raw["focus_pixel_size_unbinned_nm"],
                                "shift_x_um": raw["shift_x_um"],
                                "shift_y_um": raw["shift_y_um"],
                                "shift_abs_um": raw["shift_abs_um"],
                                "beta_rad": legacy["beta_rad"],
                                "linear_term_um": legacy["linear_term_um"],
                                "cs_term_um": legacy["cs_term_um"],
                                "spherical_aberration_mm": legacy["spherical_aberration_mm"],
                                "legacy_defocus_um": legacy_defocus,
                                "empirical_defocus_factor": float(empirical_defocus_factor),
                                "empirical_defocus_um": empirical_defocus,
                                "delta_empirical_minus_ctf_um": delta_empirical_ctf,
                                "ratio_empirical_over_ctf": ratio_empirical_ctf,
                                "fitted_empirical_defocus_factor": "",
                                "empirical_fitted_defocus_um": "",
                                "delta_empirical_fitted_minus_ctf_um": "",
                                "ratio_empirical_fitted_over_ctf": "",
                                "ctf_defocus_um": ctf_defocus,
                                "delta_legacy_minus_ctf_um": delta_legacy_ctf,
                                "ratio_legacy_over_ctf": ratio_legacy_ctf,
                                "ctf_resolution_A": ctf_res,
                                "drift_speed_x_nm_per_s": raw["drift_speed_x_nm_per_s"],
                                "drift_speed_y_nm_per_s": raw["drift_speed_y_nm_per_s"],
                                "drift_speed_abs_nm_per_s": drift_abs,
                                "elapsed_s": raw["elapsed_s"],
                                "fit_used": int(fit_used),
                            }
                            append_measurement(csv_measurements, row)
                            rows.append(row)
                            ratio_legacy_text = (
                                f"{ratio_legacy_ctf:.4f}"
                                if np.isfinite(ratio_legacy_ctf)
                                else "NaN"
                            )
                            ratio_empirical_text = (
                                f"{ratio_empirical_ctf:.4f}"
                                if np.isfinite(ratio_empirical_ctf)
                                else "NaN"
                            )
                            sem.Echo(
                                f"xtilt=({xtilt_x:.6f}, {xtilt_y:.6f}), "
                                f"bt_offset=({bt_off_x:.3f}, {bt_off_y:.3f}), "
                                f"legacy={legacy_defocus:.4f}, "
                                f"empirical(beta*f={empirical_defocus_factor:.3f})="
                                f"{empirical_defocus:.4f}, "
                                f"CTF={ctf_defocus:.4f}, "
                                f"delta_legacy={delta_legacy_ctf:.4f}, "
                                f"ratio_legacy={ratio_legacy_text}, "
                                f"delta_empirical={delta_empirical_ctf:.4f}, "
                                f"ratio_empirical={ratio_empirical_text}, "
                                f"drift={drift_abs:.3f} nm/s"
                            )
    finally:
        sem.SetBeamTilt(float(original_beam_tilt[0]), float(original_beam_tilt[1]))
        sem.SetXLensDeflector(2, float(original_xtilt[0]), float(original_xtilt[1]))
        restore_stage(original_stage)
    return rows


def feature_matrix(rows, feature_names):
    matrix = []
    for row in rows:
        raw = {
            "shift_x_um": float(row["shift_x_um"]),
            "shift_y_um": float(row["shift_y_um"]),
            "shift_abs_um": float(row["shift_abs_um"]),
            "xtilt_x": float(row["xtilt_x"]),
            "xtilt_y": float(row["xtilt_y"]),
            "beam_tilt_x0": float(row["beam_tilt_x0"]),
            "beam_tilt_y0": float(row["beam_tilt_y0"]),
            "tilt_step_x": float(row["tilt_step_x"]),
            "tilt_step_y": float(row["tilt_step_y"]),
        }
        matrix.append([btdef._feature_value(name, raw) for name in feature_names])
    return np.array(matrix, dtype=float)


def fit_defocus_model(rows):
    fit_rows = [row for row in rows if int(row["fit_used"]) == 1]
    if len(fit_rows) < 8:
        sem.Echo("WARNING: Not enough valid rows for fitted calibration.")
        return None, fit_rows

    feature_names = [
        "1",
        "shift_x_um",
        "shift_y_um",
        "xtilt_x",
        "xtilt_y",
        "beam_tilt_x0",
        "beam_tilt_y0",
        "shift_x_um*xtilt_x",
        "shift_x_um*xtilt_y",
        "shift_y_um*xtilt_x",
        "shift_y_um*xtilt_y",
    ]
    x = feature_matrix(fit_rows, feature_names)
    ctf = np.array([float(row["ctf_defocus_um"]) for row in fit_rows], dtype=float)
    legacy_base = np.array(
        [float(row["legacy_defocus_um"]) for row in fit_rows], dtype=float
    )
    y_residual = ctf - legacy_base
    coeffs, residuals, rank, singular_values = np.linalg.lstsq(x, y_residual, rcond=None)
    predicted_residual = x.dot(coeffs)
    fitted = legacy_base + predicted_residual
    errors = fitted - ctf
    rms = float(np.sqrt(np.mean(errors * errors)))

    legacy_errors = legacy_base - ctf
    legacy_rms = float(np.sqrt(np.mean(legacy_errors * legacy_errors)))
    raw_range_names = [
        "shift_x_um",
        "shift_y_um",
        "shift_abs_um",
        "xtilt_x",
        "xtilt_y",
        "beam_tilt_x0",
        "beam_tilt_y0",
        "tilt_step_x",
        "tilt_step_y",
    ]
    raw_ranges = {}
    for name in raw_range_names:
        values = np.array([float(row[name]) for row in fit_rows], dtype=float)
        raw_ranges[name] = [float(np.min(values)), float(np.max(values))]

    calibration = {
        "model": "linear_xtilt_residual",
        "created": datetime.now().isoformat(timespec="seconds"),
        "feature_names": feature_names,
        "coefficients": [float(v) for v in coeffs],
        "fit_target": "ctf_minus_legacy_physics_cs",
        "tilt_angle_mrad": float(tilt_angle_mrad),
        "beam_tilt_correction": float(beam_tilt_correction),
        "beam_tilt_axis": beam_tilt_axis,
        "spherical_aberration_mm": float(spherical_aberration_mm),
        "nominal_xtilt_x": float(nominal_xtilt_x),
        "nominal_xtilt_y": float(nominal_xtilt_y),
        "ctf_xtilt_x": float(ctf_xtilt_x),
        "ctf_xtilt_y": float(ctf_xtilt_y),
        "stage_z_offsets_um": [float(v) for v in stage_z_offsets_um],
        "use_stage_z_offsets": bool(use_stage_z_offsets),
        "stage_z_backlash_um": float(stage_z_backlash_um),
        "fit_row_count": len(fit_rows),
        "fit_rank": int(rank),
        "singular_values": [float(v) for v in singular_values],
        "rms_error_um": rms,
        "legacy_rms_error_um": legacy_rms,
        "ctf_resolution_max_A": float(ctf_resolution_max_A),
        "max_drift_nm_per_s": float(max_drift_nm_per_s),
        "raw_ranges": raw_ranges,
    }
    return calibration, fit_rows


def rows_for_empirical_fit(rows):
    fit_rows = [row for row in rows if int(row["fit_used"]) == 1]
    return fit_rows if fit_rows else rows


def empirical_rms_error_um(fit_rows, factor):
    if len(fit_rows) == 0 or not np.isfinite(factor):
        return np.nan
    errors = []
    for row in fit_rows:
        empirical = btdef.empirical_physics_cs_defocus(
            row_to_raw(row),
            tilt_angle_mrad=tilt_angle_mrad,
            cs_mm=spherical_aberration_mm,
            beam_tilt_axis=beam_tilt_axis,
            beta_factor=factor,
        )
        errors.append(empirical - float(row["ctf_defocus_um"]))
    return float(np.sqrt(np.mean(np.array(errors, dtype=float) ** 2)))


def fit_empirical_defocus_model(fit_rows):
    factors = list(empirical_defocus_factors_for_fit)
    if len(factors) < 2:
        factors = list(np.linspace(0.7, 1.3, 25))
    slope, intercept, factor_zero = fit_empirical_beta_factor_zero_crossing(
        fit_rows, factors
    )
    factor_min_rms, rms_min = fit_empirical_beta_factor_min_rms(fit_rows, factors)
    return {
        "trial_empirical_defocus_factor": float(empirical_defocus_factor),
        "trial_empirical_rms_error_um": empirical_rms_error_um(
            fit_rows, empirical_defocus_factor
        ),
        "fitted_empirical_defocus_factor": float(factor_min_rms),
        "fitted_empirical_rms_error_um": float(rms_min),
        "fitted_empirical_defocus_factor_zero_crossing": float(factor_zero),
        "zero_crossing_fit_slope": float(slope),
        "zero_crossing_fit_intercept": float(intercept),
        "beta_factor_grid": [float(v) for v in factors],
        "equation": "displacement/(2*beta*f) - Cs*(beta*f)^2",
    }


def fit_shift_sensitivity(rows):
    fit_rows = [row for row in rows if int(row["fit_used"]) == 1]
    if len(fit_rows) < 5:
        return []
    features = ["1", "ctf_defocus_um", "xtilt_x", "xtilt_y",
                "beam_tilt_x0", "beam_tilt_y0"]
    x = []
    for row in fit_rows:
        x.append([
            1.0,
            float(row["ctf_defocus_um"]),
            float(row["xtilt_x"]),
            float(row["xtilt_y"]),
            float(row["beam_tilt_x0"]),
            float(row["beam_tilt_y0"]),
        ])
    x = np.array(x, dtype=float)
    summaries = []
    for target in ("shift_x_um", "shift_y_um"):
        y = np.array([float(row[target]) for row in fit_rows], dtype=float)
        coeffs, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
        pred = x.dot(coeffs)
        err = pred - y
        summaries.append({
            "target": target,
            "feature_names": features,
            "coefficients": [float(v) for v in coeffs],
            "rank": int(rank),
            "rms_error_um": float(np.sqrt(np.mean(err * err))),
        })
    return summaries


def save_summary(path, calibration, shift_summaries, empirical_fit=None):
    fields = ["metric", "value"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        if calibration is None:
            writer.writerow({"metric": "fit_status", "value": "failed"})
        else:
            writer.writerow({"metric": "fit_status", "value": "ok"})
            writer.writerow({
                "metric": "fit_row_count",
                "value": calibration["fit_row_count"],
            })
            writer.writerow({
                "metric": "matrix_rms_error_um",
                "value": calibration["rms_error_um"],
            })
            writer.writerow({
                "metric": "legacy_rms_error_um",
                "value": calibration["legacy_rms_error_um"],
            })
        if empirical_fit is not None:
            writer.writerow({
                "metric": "trial_empirical_defocus_factor",
                "value": empirical_fit["trial_empirical_defocus_factor"],
            })
            writer.writerow({
                "metric": "trial_empirical_rms_error_um",
                "value": empirical_fit["trial_empirical_rms_error_um"],
            })
            writer.writerow({
                "metric": "fitted_empirical_defocus_factor",
                "value": empirical_fit["fitted_empirical_defocus_factor"],
            })
            writer.writerow({
                "metric": "fitted_empirical_rms_error_um",
                "value": empirical_fit["fitted_empirical_rms_error_um"],
            })
            writer.writerow({
                "metric": "fitted_empirical_defocus_factor_zero_crossing",
                "value": empirical_fit["fitted_empirical_defocus_factor_zero_crossing"],
            })
        for summary in shift_summaries:
            writer.writerow({
                "metric": f"{summary['target']}_rms_error_um",
                "value": summary["rms_error_um"],
            })


def make_plots(path, rows, calibration):
    valid_rows = [row for row in rows if int(row["fit_used"]) == 1]
    if len(valid_rows) == 0:
        return
    ctf = np.array([float(row["ctf_defocus_um"]) for row in valid_rows], dtype=float)
    legacy = np.array([float(row["legacy_defocus_um"]) for row in valid_rows], dtype=float)
    empirical = np.array(
        [float(row["empirical_defocus_um"]) for row in valid_rows], dtype=float
    )
    empirical_fitted = np.array(
        [
            float(row["empirical_fitted_defocus_um"])
            if row.get("empirical_fitted_defocus_um") not in ("", None)
            else np.nan
            for row in valid_rows
        ],
        dtype=float,
    )
    shift_x = np.array([float(row["shift_x_um"]) for row in valid_rows], dtype=float)
    shift_y = np.array([float(row["shift_y_um"]) for row in valid_rows], dtype=float)
    xtilt_x = np.array([float(row["xtilt_x"]) for row in valid_rows], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), tight_layout=True)
    axes = axes.ravel()
    axes[0].scatter(ctf, legacy - ctf)
    axes[0].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_xlabel("CTF defocus (um)")
    axes[0].set_ylabel("Legacy physics+Cs - CTF (um)")
    axes[0].set_title("Known physics error (before matrix correction)")

    axes[1].scatter(xtilt_x, legacy - ctf, label="legacy", alpha=0.7)
    if np.any(np.isfinite(empirical)):
        axes[1].scatter(xtilt_x, empirical - ctf, label="empirical", alpha=0.7)
    axes[1].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_xlabel("XLensDeflector(2) X")
    axes[1].set_ylabel("Defocus - CTF (um)")
    axes[1].set_title("X-tilt sensitivity")
    axes[1].legend()

    axes[2].scatter(ctf, shift_x, label="shift_x")
    axes[2].scatter(ctf, shift_y, label="shift_y")
    axes[2].set_xlabel("CTF defocus (um)")
    axes[2].set_ylabel("Raw shift (um)")
    axes[2].set_title("Raw shift vs true defocus")
    axes[2].legend()

    if calibration is not None:
        x = feature_matrix(valid_rows, calibration["feature_names"])
        residual = x.dot(np.array(calibration["coefficients"], dtype=float))
        fitted = legacy + residual
        axes[3].scatter(ctf, fitted - ctf, label="matrix", alpha=0.7)
        if np.any(np.isfinite(empirical_fitted)):
            axes[3].scatter(
                ctf, empirical_fitted - ctf, label="empirical fitted", alpha=0.7
            )
        axes[3].axhline(0, color="gray", linestyle="--", linewidth=1)
        axes[3].set_xlabel("CTF defocus (um)")
        axes[3].set_ylabel("Fitted - CTF (um)")
        axes[3].set_title("Calibrated error")
        axes[3].legend()
    else:
        axes[3].axis("off")

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    sem.Echo("##### Beam tilt / X-tilt matrix calibration #####")
    sem.Echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    prepare_output_paths()
    if save_dir:
        sem.Echo(f"Output directory: {save_dir}")
    initialize_measurements(csv_measurements)
    rows = collect_measurements()
    calibration, fit_rows = fit_defocus_model(rows)
    shift_summaries = fit_shift_sensitivity(rows)
    empirical_fit_rows = rows_for_empirical_fit(rows)
    empirical_fit = fit_empirical_defocus_model(empirical_fit_rows)
    fitted_empirical_factor = empirical_fit["fitted_empirical_defocus_factor"]
    if np.isfinite(fitted_empirical_factor):
        apply_empirical_fitted_columns(rows, fitted_empirical_factor)
        write_measurements_csv(csv_measurements, rows)

    if calibration is not None:
        calibration["shift_sensitivity_models"] = shift_summaries
        calibration["empirical_beta_factor_fit"] = empirical_fit
        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)
        sem.Echo(
            f"Matrix RMS error: {calibration['rms_error_um']:.4f} um "
            f"(legacy physics+Cs {calibration['legacy_rms_error_um']:.4f} um)"
        )
        sem.Echo(f"Saved calibration JSON: {calibration_json}")
    else:
        with open(calibration_json, "w") as fh:
            json.dump({"empirical_beta_factor_fit": empirical_fit}, fh, indent=2)
        sem.Echo(f"Saved empirical-only JSON: {calibration_json}")

    sem.Echo(
        f"Empirical trial factor {empirical_fit['trial_empirical_defocus_factor']:.4f}: "
        f"RMS {empirical_fit['trial_empirical_rms_error_um']:.4f} um"
    )
    sem.Echo(
        f"Empirical beta factor (min RMS) "
        f"{empirical_fit['fitted_empirical_defocus_factor']:.4f}: "
        f"RMS {empirical_fit['fitted_empirical_rms_error_um']:.4f} um"
    )
    sem.Echo(
        f"Empirical fitted factor (zero crossing) "
        f"{empirical_fit['fitted_empirical_defocus_factor_zero_crossing']:.4f}"
    )

    save_summary(csv_summary, calibration, shift_summaries, empirical_fit)
    make_plots(plot_file, rows, calibration)
    sem.Echo(f"Saved measurements CSV: {csv_measurements}")
    sem.Echo(f"Saved summary CSV: {csv_summary}")
    sem.Echo(f"Saved plot: {plot_file}")
    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
