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
tilt_angle_mrad = 5.0
beam_tilt_correction = 1.73
beam_tilt_axis = "x"

# Nominal X-tilt used as the center of the beam-tilt autofocus sweep.
nominal_xtilt_x = 0.0
nominal_xtilt_y = 0.0

# X-tilt used for CTF reference defocus. This can be the high/known-good
# PACEtomo.py ctfXtiltX/Y setting, independent of the X-tilt being tested.
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867

# Test X-tilt around nominal. Start one axis at a time; enable the 2D grid
# after the one-axis tests show the sweep range is safe.
xtilt_step = 0.00005
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
target_defocus_tolerance_um = 0.02
max_defocus_adjust_iterations = 5

# Repeats per condition.
repeats = 1

# Reject rows from fitting if drift exceeds this speed. Set <= 0 to disable.
max_drift_nm_per_s = 2.0

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = ""

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
        "shift_x_um",
        "shift_y_um",
        "shift_abs_um",
        "legacy_defocus_um",
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
                            legacy_defocus = btdef.defocus_from_raw(
                                raw,
                                tilt_angle_mrad=tilt_angle_mrad,
                                legacy_divisor=2.0,
                                calibration_data={"model": "legacy_radial"},
                            )
                            ctf_defocus, ctf_res = acquire_ctf_reference()
                            delta_legacy_ctf = legacy_defocus - ctf_defocus
                            ratio_legacy_ctf = legacy_ctf_ratio(
                                legacy_defocus, ctf_defocus
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
                                "shift_x_um": raw["shift_x_um"],
                                "shift_y_um": raw["shift_y_um"],
                                "shift_abs_um": raw["shift_abs_um"],
                                "legacy_defocus_um": legacy_defocus,
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
                            ratio_text = (
                                f"{ratio_legacy_ctf:.4f}"
                                if np.isfinite(ratio_legacy_ctf)
                                else "NaN"
                            )
                            sem.Echo(
                                f"xtilt=({xtilt_x:.6f}, {xtilt_y:.6f}), "
                                f"bt_offset=({bt_off_x:.3f}, {bt_off_y:.3f}), "
                                f"legacy={legacy_defocus:.4f}, CTF={ctf_defocus:.4f}, "
                                f"delta={delta_legacy_ctf:.4f}, "
                                f"ratio={ratio_text}, "
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
    y = np.array([float(row["ctf_defocus_um"]) for row in fit_rows], dtype=float)
    coeffs, residuals, rank, singular_values = np.linalg.lstsq(x, y, rcond=None)
    predicted = x.dot(coeffs)
    errors = predicted - y
    rms = float(np.sqrt(np.mean(errors * errors)))

    legacy_errors = np.array(
        [float(row["legacy_defocus_um"]) - float(row["ctf_defocus_um"])
         for row in fit_rows],
        dtype=float,
    )
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
        "model": "linear_xtilt_beam_tilt",
        "created": datetime.now().isoformat(timespec="seconds"),
        "feature_names": feature_names,
        "coefficients": [float(v) for v in coeffs],
        "tilt_angle_mrad": float(tilt_angle_mrad),
        "beam_tilt_correction": float(beam_tilt_correction),
        "beam_tilt_axis": beam_tilt_axis,
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


def save_summary(path, calibration, shift_summaries):
    fields = ["metric", "value"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        if calibration is None:
            writer.writerow({"metric": "fit_status", "value": "failed"})
            return
        writer.writerow({"metric": "fit_status", "value": "ok"})
        writer.writerow({"metric": "fit_row_count", "value": calibration["fit_row_count"]})
        writer.writerow({"metric": "rms_error_um", "value": calibration["rms_error_um"]})
        writer.writerow({
            "metric": "legacy_rms_error_um",
            "value": calibration["legacy_rms_error_um"],
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
    shift_x = np.array([float(row["shift_x_um"]) for row in valid_rows], dtype=float)
    shift_y = np.array([float(row["shift_y_um"]) for row in valid_rows], dtype=float)
    xtilt_x = np.array([float(row["xtilt_x"]) for row in valid_rows], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), tight_layout=True)
    axes = axes.ravel()
    axes[0].scatter(ctf, legacy - ctf)
    axes[0].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_xlabel("CTF defocus (um)")
    axes[0].set_ylabel("Legacy - CTF (um)")
    axes[0].set_title("Legacy beam-tilt error")

    axes[1].scatter(xtilt_x, legacy - ctf)
    axes[1].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_xlabel("XLensDeflector(2) X")
    axes[1].set_ylabel("Legacy - CTF (um)")
    axes[1].set_title("X-tilt sensitivity")

    axes[2].scatter(ctf, shift_x, label="shift_x")
    axes[2].scatter(ctf, shift_y, label="shift_y")
    axes[2].set_xlabel("CTF defocus (um)")
    axes[2].set_ylabel("Raw shift (um)")
    axes[2].set_title("Raw shift vs true defocus")
    axes[2].legend()

    if calibration is not None:
        x = feature_matrix(valid_rows, calibration["feature_names"])
        fitted = x.dot(np.array(calibration["coefficients"], dtype=float))
        axes[3].scatter(ctf, fitted - ctf)
        axes[3].axhline(0, color="gray", linestyle="--", linewidth=1)
        axes[3].set_xlabel("CTF defocus (um)")
        axes[3].set_ylabel("Fitted - CTF (um)")
        axes[3].set_title("Fitted calibration residual")
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

    if calibration is not None:
        calibration["shift_sensitivity_models"] = shift_summaries
        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)
        sem.Echo(
            f"Fitted RMS error: {calibration['rms_error_um']:.4f} um "
            f"(legacy {calibration['legacy_rms_error_um']:.4f} um)"
        )
        sem.Echo(f"Saved calibration JSON: {calibration_json}")

    save_summary(csv_summary, calibration, shift_summaries)
    make_plots(plot_file, rows, calibration)
    sem.Echo(f"Saved measurements CSV: {csv_measurements}")
    sem.Echo(f"Saved summary CSV: {csv_summary}")
    sem.Echo(f"Saved plot: {plot_file}")
    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
