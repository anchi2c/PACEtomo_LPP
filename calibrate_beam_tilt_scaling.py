#!Python
# ===================================================================
# ScriptName     calibrate_beam_tilt_scaling
# Purpose:       At X-tilt 0, sweep defocus and fit delta/ratio vs CTF to
#                compensate beam-tilt defocus (physics + Cs).
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

tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
beam_tilt_axis = "x"

target_defocus_values = [-1.0, -2.0, -3.0, -4.0, -5.0]

ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2
target_defocus_tolerance_um = 0.05
max_defocus_adjust_iterations = 5

xtilt_lens_index = 2
measurement_xtilt_x = 0.0
measurement_xtilt_y = 0.0

save_dir = ""
csv_measurements = "beam_tilt_scaling_measurements.csv"
calibration_json = "beam_tilt_scaling_calibration.json"
plot_file = "beam_tilt_scaling_fits.png"
shift_plot_file = "beam_tilt_scaling_shifts.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


btdef.configure(sem_module=sem, logger=echo)


def prepare_output_paths():
    global csv_measurements, calibration_json, plot_file, shift_plot_file
    if not save_dir:
        return
    os.makedirs(save_dir, exist_ok=True)
    csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
    calibration_json = os.path.join(save_dir, os.path.basename(calibration_json))
    plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    shift_plot_file = os.path.join(save_dir, os.path.basename(shift_plot_file))


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


def set_target_defocus(target_defocus):
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    current_defocus = np.nan
    for attempt in range(1, max_defocus_adjust_iterations + 1):
        sem.L()
        current_defocus, _ = run_ctffind()
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


def measure_beam_tilt_defocus():
    raw = btdef.measure_raw(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        beam_tilt_axis=beam_tilt_axis,
    )
    diag = btdef.legacy_physics_diagnostics(
        raw,
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_axis=beam_tilt_axis,
        defocus_tilt_correction=beam_tilt_correction,
    )
    return {
        "beam_tilt_defocus_um": float(diag["legacy_defocus_um"]),
        "cs_term_um": float(diag["cs_term_um"]),
        "beta_rad": float(diag["beta_rad"]),
        "raw": raw,
    }


def safe_ratio(measured, ctf):
    if not np.isfinite(measured) or not np.isfinite(ctf) or abs(ctf) < 1e-6:
        return np.nan
    return float(measured) / float(ctf)


def fit_line(x, y):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.sum(mask) < 2:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    pred = slope * x[mask] + intercept
    rms = float(np.sqrt(np.mean((pred - y[mask]) ** 2)))
    return float(slope), float(intercept), rms


def fit_affine_true_from_measured(measured, ctf):
    measured = np.array(measured, dtype=float)
    ctf = np.array(ctf, dtype=float)
    mask = np.isfinite(measured) & np.isfinite(ctf)
    if np.sum(mask) < 2:
        return np.nan, np.nan, np.nan
    x = measured[mask]
    y = ctf[mask]
    scale, offset = np.polyfit(x, y, 1)
    pred = scale * x + offset
    rms = float(np.sqrt(np.mean((pred - y) ** 2)))
    return float(scale), float(offset), rms


def fit_scaling(measured, ctf):
    measured = [float(v) for v in measured]
    ctf = [float(v) for v in ctf]
    deltas = [m - c for m, c in zip(measured, ctf)]
    ratios = [safe_ratio(m, c) for m, c in zip(measured, ctf)]

    delta_slope, delta_intercept, delta_rms = fit_line(ctf, deltas)
    ratio_slope, ratio_intercept, ratio_rms = fit_line(ctf, ratios)
    affine_scale, affine_offset, affine_rms = fit_affine_true_from_measured(
        measured, ctf
    )
    return {
        "delta": {
            "intercept": delta_intercept,
            "slope_vs_ctf": delta_slope,
            "rms_um": delta_rms,
            "formula": "delta = intercept + slope * ctf_defocus",
        },
        "ratio": {
            "intercept": ratio_intercept,
            "slope_vs_ctf": ratio_slope,
            "rms": ratio_rms,
            "formula": "ratio = intercept + slope * ctf_defocus",
        },
        "affine": {
            "scale": affine_scale,
            "offset": affine_offset,
            "rms_um": affine_rms,
            "formula": "corrected_defocus = scale * measured + offset",
        },
        "constant": {
            "mean_delta_um": float(np.nanmean(deltas)),
            "mean_ratio": float(np.nanmean(ratios)),
        },
    }


def make_scaling_plot(path, rows, fits):
    ctf = np.array([float(r["ctf_defocus_um"]) for r in rows], dtype=float)
    measured = np.array([float(r["beam_tilt_defocus_um"]) for r in rows], dtype=float)
    delta = np.array([float(r["delta_um"]) for r in rows], dtype=float)
    ratio = np.array([float(r["ratio"]) for r in rows], dtype=float)
    ctf_line = np.linspace(np.min(ctf), np.max(ctf), 100)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), tight_layout=True)

    ax = axes[0, 0]
    ax.scatter(ctf, measured, label="data")
    scale = fits["affine"]["scale"]
    offset = fits["affine"]["offset"]
    if np.isfinite(scale) and scale != 0:
        ax.plot(ctf_line, (ctf_line - offset) / scale, label="affine", color="C1")
    ax.plot(ctf_line, ctf_line, "--", color="gray", label="ideal")
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Beam-tilt defocus (um)")
    ax.set_title("Measured vs CTF")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.scatter(ctf, delta)
    ax.plot(
        ctf_line,
        fits["delta"]["intercept"] + fits["delta"]["slope_vs_ctf"] * ctf_line,
        color="C1",
    )
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Delta (um)")
    ax.set_title("Delta vs CTF")

    ax = axes[1, 0]
    ax.scatter(ctf, ratio)
    ax.plot(
        ctf_line,
        fits["ratio"]["intercept"] + fits["ratio"]["slope_vs_ctf"] * ctf_line,
        color="C1",
    )
    ax.axhline(1, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Ratio")
    ax.set_title("Ratio vs CTF")

    corrected = np.array(
        [
            fits["affine"]["scale"] * r["beam_tilt_defocus_um"]
            + fits["affine"]["offset"]
            for r in rows
        ],
        dtype=float,
    )
    ax = axes[1, 1]
    ax.scatter(ctf, corrected - ctf)
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Corrected - CTF (um)")
    ax.set_title("Affine correction residual")

    fig.savefig(path, dpi=150)
    plt.show()


def make_shift_plot(path, rows, plus_fit, minus_fit):
    ctf = np.array([float(r["ctf_defocus_um"]) for r in rows], dtype=float)
    plus = np.array([float(r["plus_branch_shift_um"]) for r in rows], dtype=float)
    minus = np.array([float(r["minus_branch_shift_um"]) for r in rows], dtype=float)
    ctf_line = np.linspace(np.min(ctf), np.max(ctf), 100)

    fig, ax = plt.subplots(figsize=(8, 5), tight_layout=True)
    ax.scatter(ctf, plus, label="+ branch", color="C0")
    ax.scatter(ctf, minus, label="- branch", color="C1")
    ax.plot(
        ctf_line,
        plus_fit[1] + plus_fit[0] * ctf_line,
        color="C0",
        label=f"+ fit: {plus_fit[1]:.4f} + {plus_fit[0]:.4f}*ctf",
    )
    ax.plot(
        ctf_line,
        minus_fit[1] + minus_fit[0] * ctf_line,
        color="C1",
        label=f"- fit: {minus_fit[1]:.4f} + {minus_fit[0]:.4f}*ctf",
    )
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel(f"Signed shift along {beam_tilt_axis} (um)")
    ax.set_title("Plus/minus branch shifts vs defocus")
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=150)
    plt.show()


def echo_fit_summary(fits):
    echo(f"  mean delta: {fits['constant']['mean_delta_um']:.6f} um")
    echo(f"  mean ratio: {fits['constant']['mean_ratio']:.6f}")
    echo(
        f"  delta vs CTF: {fits['delta']['intercept']:.6f} "
        f"+ {fits['delta']['slope_vs_ctf']:.6f} * ctf  "
        f"(RMS {fits['delta']['rms_um']:.6f} um)"
    )
    echo(
        f"  ratio vs CTF: {fits['ratio']['intercept']:.6f} "
        f"+ {fits['ratio']['slope_vs_ctf']:.6f} * ctf  "
        f"(RMS {fits['ratio']['rms']:.6f})"
    )
    echo(
        f"  affine: corrected = {fits['affine']['scale']:.6f} * measured "
        f"+ {fits['affine']['offset']:.6f}  "
        f"(RMS {fits['affine']['rms_um']:.6f} um)"
    )


def main():
    sem.SuppressReports()
    prepare_output_paths()
    echo("##### Beam tilt defocus scaling calibration (X-tilt 0) #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    if save_dir:
        echo(f"Output directory: {save_dir}")
    echo(
        f"beam_tilt_correction={beam_tilt_correction}, "
        f"tilt_angle_mrad={tilt_angle_mrad}, "
        f"SetBeamTilt step={beam_tilt_correction * tilt_angle_mrad:.4f}"
    )

    original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    echo(
        f"Saved XLensDeflector({xtilt_lens_index}): "
        f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
    )

    rows = []
    try:
        sem.SetXLensDeflector(
            xtilt_lens_index, float(measurement_xtilt_x), float(measurement_xtilt_y)
        )
        echo(
            f"Set XLensDeflector({xtilt_lens_index}) to "
            f"({measurement_xtilt_x:.6f}, {measurement_xtilt_y:.6f})"
        )

        for target_defocus in target_defocus_values:
            echo("------------------------------------------------")
            echo(f"Target defocus {float(target_defocus):.3f} um")
            reached = set_target_defocus(target_defocus)
            echo(f"Defocus after adjustment: {reached:.3f} um")

            measure = measure_beam_tilt_defocus()
            raw = measure["raw"]
            sem.L()
            ctf_defocus, ctf_res = run_ctffind()
            drift = float(np.hypot(
                raw["drift_speed_x_nm_per_s"], raw["drift_speed_y_nm_per_s"]
            ))

            row = {
                "target_defocus_um": float(target_defocus),
                "ctf_defocus_um": float(ctf_defocus),
                "beam_tilt_defocus_um": measure["beam_tilt_defocus_um"],
                "cs_term_um": measure["cs_term_um"],
                "beta_rad": measure["beta_rad"],
                "delta_um": measure["beam_tilt_defocus_um"] - ctf_defocus,
                "ratio": safe_ratio(measure["beam_tilt_defocus_um"], ctf_defocus),
                "shift_axis_um": float(raw["shift_axis_um"]),
                "plus_branch_shift_um": float(raw["plus_branch_shift_um"]),
                "minus_branch_shift_um": float(raw["minus_branch_shift_um"]),
                "minus_vs_plus_ref_um": float(raw["minus_vs_plus_ref_um"]),
                "drift_return_um": float(raw["drift_return_um"]),
                "shift_x_um": float(raw["shift_x_um"]),
                "shift_y_um": float(raw["shift_y_um"]),
                "ctf_resolution_A": float(ctf_res),
                "drift_speed_nm_per_s": drift,
            }
            rows.append(row)
            echo(
                f"defocus={measure['beam_tilt_defocus_um']:.4f}, "
                f"Cs term={measure['cs_term_um']:.4f}, "
                f"CTF={ctf_defocus:.4f} um, "
                f"delta={row['delta_um']:.4f}, ratio={row['ratio']:.4f}, "
                f"+ branch={row['plus_branch_shift_um']:.4f} um, "
                f"- branch={row['minus_branch_shift_um']:.4f} um, "
                f"- vs + ref={row['minus_vs_plus_ref_um']:.4f} um, "
                f"drift={drift:.3f} nm/s"
            )

        ctf = [r["ctf_defocus_um"] for r in rows]
        plus = [r["plus_branch_shift_um"] for r in rows]
        minus = [r["minus_branch_shift_um"] for r in rows]
        fits = fit_scaling([r["beam_tilt_defocus_um"] for r in rows], ctf)
        plus_slope, plus_intercept, plus_rms = fit_line(ctf, plus)
        minus_slope, minus_intercept, minus_rms = fit_line(ctf, minus)

        calibration = {
            "model": "beam_tilt_defocus_scaling",
            "created": datetime.now().isoformat(timespec="seconds"),
            "tilt_angle_mrad": float(tilt_angle_mrad),
            "beam_tilt_correction": float(beam_tilt_correction),
            "beam_tilt_axis": beam_tilt_axis,
            "spherical_aberration_mm": float(btdef.spherical_aberration_mm),
            "measurement_xtilt_x": float(measurement_xtilt_x),
            "measurement_xtilt_y": float(measurement_xtilt_y),
            "equation": "-displacement/(2*beta) - Cs*beta^2",
            "fits": fits,
            "branch_shift_fits_vs_ctf": {
                "plus_branch": {
                    "slope": plus_slope,
                    "intercept": plus_intercept,
                    "rms_um": plus_rms,
                    "formula": "plus_shift = intercept + slope * ctf_defocus",
                },
                "minus_branch": {
                    "slope": minus_slope,
                    "intercept": minus_intercept,
                    "rms_um": minus_rms,
                    "formula": "minus_shift = intercept + slope * ctf_defocus",
                },
            },
            "recommended": {
                "model": "affine",
                "scale": fits["affine"]["scale"],
                "offset": fits["affine"]["offset"],
                "rms_um": fits["affine"]["rms_um"],
            },
        }

        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)

        make_scaling_plot(plot_file, rows, fits)
        make_shift_plot(
            shift_plot_file, rows, (plus_slope, plus_intercept), (minus_slope, minus_intercept)
        )

        echo("================================================")
        echo("FIT RESULTS (X-tilt 0, physics + Cs)")
        echo_fit_summary(fits)
        echo(
            f"+ branch vs CTF: {plus_intercept:.6f} + {plus_slope:.6f} * ctf "
            f"(RMS {plus_rms:.6f} um)"
        )
        echo(
            f"- branch vs CTF: {minus_intercept:.6f} + {minus_slope:.6f} * ctf "
            f"(RMS {minus_rms:.6f} um)"
        )
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved scaling plot: {plot_file}")
        echo(f"Saved shift plot: {shift_plot_file}")

    finally:
        sem.SetXLensDeflector(
            xtilt_lens_index, float(original_xtilt[0]), float(original_xtilt[1])
        )
        echo(
            f"Restored XLensDeflector({xtilt_lens_index}) to "
            f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
        )

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
