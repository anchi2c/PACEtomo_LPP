#!Python
# ===================================================================
# ScriptName     calibrate_beam_tilt_scaling
# Purpose:       Sweep defocus at a fixed X-tilt and fit the systematic
#                offset between physics+Cs beam-tilt defocus and CTF:
#                  delta = measured - ctf = intercept + slope * defocus
#                Writes beam_tilt_scaling_calibration.json for
#                PACEtomo_beamTiltDefocus.calibration_file.
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

# XLensDeflector index and X-tilt for beam-tilt measurements (restored after run).
xtilt_lens_index = 2
measurement_xtilt_x = 0.0
measurement_xtilt_y = 0.0

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r"X:\k3f_leginonframes\p26jun29a\xtilt_calib_test"

# Output file names (written under save_dir when set).
csv_measurements = "beam_tilt_scaling_measurements.csv"
calibration_json = "beam_tilt_scaling_calibration.json"
plot_file = "beam_tilt_scaling_fits.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


btdef.configure(sem_module=sem, logger=echo)


def prepare_output_paths():
    """Resolve output paths under save_dir and create the directory if needed."""
    global csv_measurements, calibration_json, plot_file
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
        calibration_json = os.path.join(save_dir, os.path.basename(calibration_json))
        plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    csv_measurements = os.path.abspath(csv_measurements)
    calibration_json = os.path.abspath(calibration_json)
    plot_file = os.path.abspath(plot_file)


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
        "drift_speed_nm_per_s": float(np.hypot(
            raw["drift_speed_x_nm_per_s"], raw["drift_speed_y_nm_per_s"]
        )),
    }


def fit_delta_vs_defocus(defocus_um, delta_um):
    x = np.array(defocus_um, dtype=float)
    y = np.array(delta_um, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.sum(mask) < 2:
        return {
            "intercept_um": np.nan,
            "slope": np.nan,
            "rms_um": np.nan,
            "mean_delta_um": float(np.nanmean(y)),
        }
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    pred = slope * x[mask] + intercept
    return {
        "intercept_um": float(intercept),
        "slope": float(slope),
        "rms_um": float(np.sqrt(np.mean((pred - y[mask]) ** 2))),
        "mean_delta_um": float(np.nanmean(y[mask])),
    }


def make_plot(path, rows, delta_fit):
    ctf = np.array([float(r["ctf_defocus_um"]) for r in rows], dtype=float)
    measured = np.array([float(r["beam_tilt_defocus_um"]) for r in rows], dtype=float)
    delta = np.array([float(r["delta_um"]) for r in rows], dtype=float)
    ctf_line = np.linspace(np.min(ctf), np.max(ctf), 100)
    delta_line = delta_fit["intercept_um"] + delta_fit["slope"] * ctf_line
    corrected = measured - (delta_fit["intercept_um"] + delta_fit["slope"] * ctf)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)

    ax = axes[0]
    ax.scatter(ctf, measured, label="measured")
    ax.plot(ctf_line, ctf_line, "--", color="gray", label="ideal")
    ax.scatter(ctf, corrected, marker="x", label="corrected")
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Beam-tilt defocus (um)")
    ax.set_title("Measured and corrected vs CTF")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(ctf, delta, label="data")
    ax.plot(ctf_line, delta_line, color="C1", label="fit")
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CTF defocus (um)")
    ax.set_ylabel("Delta = measured - CTF (um)")
    ax.set_title("Delta offset vs defocus")
    ax.legend(fontsize=8)

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    xtilt_label = (
        f"X-tilt ({measurement_xtilt_x:.6f}, {measurement_xtilt_y:.6f})"
    )
    echo(f"##### Beam tilt defocus scaling calibration ({xtilt_label}) #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"beam_tilt_correction={beam_tilt_correction}, "
        f"tilt_angle_mrad={tilt_angle_mrad}, "
        f"SetBeamTilt step={beam_tilt_correction * tilt_angle_mrad:.4f}, "
        f"measurement {xtilt_label}"
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
            sem.L()
            ctf_defocus, ctf_res = run_ctffind()

            row = {
                "target_defocus_um": float(target_defocus),
                "ctf_defocus_um": float(ctf_defocus),
                "beam_tilt_defocus_um": measure["beam_tilt_defocus_um"],
                "cs_term_um": measure["cs_term_um"],
                "beta_rad": measure["beta_rad"],
                "delta_um": measure["beam_tilt_defocus_um"] - ctf_defocus,
                "ctf_resolution_A": float(ctf_res),
                "drift_speed_nm_per_s": measure["drift_speed_nm_per_s"],
            }
            rows.append(row)
            echo(
                f"measured={measure['beam_tilt_defocus_um']:.4f} um, "
                f"Cs term={measure['cs_term_um']:.4f} um, "
                f"CTF={ctf_defocus:.4f} um, "
                f"delta={row['delta_um']:.4f} um, "
                f"drift={measure['drift_speed_nm_per_s']:.3f} nm/s"
            )

        ctf = [r["ctf_defocus_um"] for r in rows]
        deltas = [r["delta_um"] for r in rows]
        delta_fit = fit_delta_vs_defocus(ctf, deltas)

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
            "delta_offset": {
                "intercept_um": delta_fit["intercept_um"],
                "slope": delta_fit["slope"],
                "rms_um": delta_fit["rms_um"],
                "mean_delta_um": delta_fit["mean_delta_um"],
                "formula": (
                    "correction_um = intercept_um + slope * defocus_um; "
                    "corrected_defocus_um = measured_defocus_um - correction_um"
                ),
            },
        }

        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)

        make_plot(plot_file, rows, delta_fit)

        echo("================================================")
        echo(f"CALIBRATION ({xtilt_label})")
        echo(f"  mean delta: {delta_fit['mean_delta_um']:.6f} um")
        echo(
            f"  delta vs defocus: {delta_fit['intercept_um']:.6f} "
            f"+ {delta_fit['slope']:.6f} * defocus  "
            f"(RMS {delta_fit['rms_um']:.6f} um)"
        )
        echo(
            "  Apply in PACEtomo_beamTiltDefocus: set calibration_file to "
            f"{calibration_json}"
        )
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved plot: {plot_file}")

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
