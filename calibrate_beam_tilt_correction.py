#!Python
# ===================================================================
# ScriptName     calibrate_beam_tilt_correction
# Purpose:       Calibrate beam tilt correction factor by comparing
#                autofocus-derived defocus with CtfFind defocus.
# Author:        Generated for PACEtomo workflow
# ===================================================================

import csv
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

import PACEtomo_beamTiltDefocus as btdef

############ SETTINGS ############

# Beam tilt used for autofocus measurement [mrad]
tilt_angle_mrad = 5.0

# Correction factors to test. Example: np.linspace(0.2, 0.8, 13)
beam_tilt_corrections = np.linspace(0.4, 0.5, 11)

# Target defocus values [microns] to calibrate at
target_defocus_values = [-1.0, -2.0, -3.0, -4.0, -5.0]

# CTF search range [microns]
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

# Defocus convergence settings for setting each target defocus
target_defocus_tolerance_um = 0.01
max_defocus_adjust_iterations = 5

# Number of autofocus cycles per correction test.
# Delta is computed using the last cycle's autofocus value.
autofocus_cycles = 2

# Output file names (written in current SerialEM working directory)
csv_measurements = "beam_tilt_correction_measurements.csv"
csv_summary = "beam_tilt_correction_summary.csv"
plot_file = "beam_tilt_correction_fits.png"

########## END SETTINGS ##########


def run_ctffind():
    """Run CtfFind on buffer A and return defocus in microns."""
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
    """Adjust focus to target defocus and verify with CtfFind."""
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    current_defocus = np.nan
    for attempt in range(1, max_defocus_adjust_iterations + 1):
        sem.L()
        current_defocus, _ = run_ctffind()
        error = target_defocus - current_defocus
        sem.Echo(
            f"Target defocus check {attempt}/{max_defocus_adjust_iterations}: "
            f"current={current_defocus:.3f} um, target={target_defocus:.3f} um, "
            f"error={error:.3f} um"
        )
        if abs(error) <= target_defocus_tolerance_um:
            sem.Echo("Target defocus reached within tolerance.")
            return current_defocus
        sem.ChangeFocus(error)
    sem.Echo(
        "WARNING: Target defocus not reached within tolerance after max iterations."
    )
    return current_defocus


def run_autofocus_trial(correction):
    """
    Execute one autofocus-like measurement with a given beam tilt correction.
    Returns:
      defocus_measured [microns], speed_x [nm/s], speed_y [nm/s]
    """
    beam_tilt = sem.ReportBeamTilt()
    tilt_x_orig = float(beam_tilt[0])
    tilt_y_orig = float(beam_tilt[1])
    tilt_x_plus = tilt_x_orig + correction * tilt_angle_mrad
    tilt_x_minus = tilt_x_orig - correction * tilt_angle_mrad

    focus_camera = "F"
    pixel_size_binned_nm = float(sem.ReportCurrentPixelSize(focus_camera))
    focus_binning = float(sem.ReportBinning(focus_camera))
    pixel_size_unbinned_nm = pixel_size_binned_nm / focus_binning

    # Positive tilt
    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    sem.F()
    sem.ResetClock()
    sem.Copy("A", "L")

    # Negative tilt
    sem.SetBeamTilt(tilt_x_minus, tilt_y_orig)
    sem.F()
    sem.AlignTo("L", 1)
    align_shift_1 = sem.ReportAlignShift()
    disp_x1_px = float(align_shift_1[0])
    disp_y1_px = float(align_shift_1[1])

    # Positive tilt again
    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    sem.F()
    elapsed = float(sem.ReportClock())

    # Back to origin
    sem.SetBeamTilt(tilt_x_orig, tilt_y_orig)
    sem.AlignTo("L", 1)
    align_shift_2 = sem.ReportAlignShift()
    disp_x2_px = float(align_shift_2[0])
    disp_y2_px = float(align_shift_2[1])

    drift_x_nm = disp_x2_px * pixel_size_unbinned_nm
    drift_y_nm = disp_y2_px * pixel_size_unbinned_nm
    speed_x = drift_x_nm / elapsed if elapsed > 0 else 0.0
    speed_y = drift_y_nm / elapsed if elapsed > 0 else 0.0

    shift_x_um = (disp_x1_px - disp_x2_px / 2.0) * pixel_size_unbinned_nm / 1000.0
    shift_y_um = (disp_y1_px - disp_y2_px / 2.0) * pixel_size_unbinned_nm / 1000.0
    raw = {
        "shift_x_um": shift_x_um,
        "shift_y_um": shift_y_um,
        "shift_abs_um": float(np.sqrt(shift_x_um * shift_x_um + shift_y_um * shift_y_um)),
    }
    defocus_measured = btdef.legacy_physics_diagnostics(
        raw, tilt_angle_mrad=tilt_angle_mrad
    )["legacy_defocus_um"]

    return defocus_measured, speed_x, speed_y


def fit_best_correction(corrections, deltas):
    """Linear fit of delta vs correction; returns slope, intercept, root."""
    slope, intercept = np.polyfit(corrections, deltas, 1)
    if slope == 0:
        return slope, intercept, np.nan
    return slope, intercept, -intercept / slope


def save_measurements_csv(path, rows):
    fields = [
        "target_defocus_um",
        "beam_tilt_correction",
        "autofocus_cycle_used",
        "autofocus_defocus_cycle1_um",
        "autofocus_defocus_cycle2_um",
        "autofocus_defocus_um",
        "ctf_defocus_um",
        "delta_um",
        "ctf_resolution_A",
        "drift_speed_x_nm_per_s",
        "drift_speed_y_nm_per_s",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def save_summary_csv(path, rows, mean_corr):
    fields = [
        "target_defocus_um",
        "fitted_correction",
        "fit_slope",
        "fit_intercept",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow(
            {
                "target_defocus_um": "MEAN",
                "fitted_correction": mean_corr,
                "fit_slope": "",
                "fit_intercept": "",
            }
        )


def make_plot(path, grouped_results):
    n = len(grouped_results)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), tight_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for i, res in enumerate(grouped_results):
        ax = axes[i]
        x = np.array(res["corrections"], dtype=float)
        y = np.array(res["deltas"], dtype=float)
        slope = res["slope"]
        intercept = res["intercept"]
        best = res["best_correction"]

        ax.scatter(x, y, label="measurements")
        xfit = np.linspace(np.min(x), np.max(x), 200)
        ax.plot(xfit, slope * xfit + intercept, label="linear fit")
        ax.axhline(0, color="gray", linestyle="--", linewidth=1)
        if np.isfinite(best):
            ax.axvline(best, color="red", linestyle=":", linewidth=1)
            title = f"Defocus {res['target_defocus']:.2f} um, best={best:.4f}"
        else:
            title = f"Defocus {res['target_defocus']:.2f} um, best=NaN"
        ax.set_title(title)
        ax.set_xlabel("Beam tilt correction")
        ax.set_ylabel("Autofocus - CTF defocus (um)")
        ax.legend()

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    sem.Echo("##### Beam tilt correction calibration #####")
    sem.Echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    sem.Echo(f"Testing corrections: {list(np.array(beam_tilt_corrections, dtype=float))}")
    sem.Echo(f"Target defocus values: {target_defocus_values}")

    measurement_rows = []
    summary_rows = []
    grouped_results = []
    fitted_corrections = []

    for target_defocus in target_defocus_values:
        sem.Echo("------------------------------------------------")
        sem.Echo(f"Starting target defocus {target_defocus:.3f} um")
        reached_defocus = set_target_defocus(target_defocus)
        sem.Echo(f"Defocus after adjustment loop: {reached_defocus:.3f} um")

        corr_vals = []
        delta_vals = []
        for corr in beam_tilt_corrections:
            corr = float(corr)
            cycle_defocus = []
            speed_x = 0.0
            speed_y = 0.0
            for _ in range(autofocus_cycles):
                af_defocus, speed_x, speed_y = run_autofocus_trial(corr)
                cycle_defocus.append(af_defocus)
            af_defocus_used = cycle_defocus[-1]
            cycle1 = cycle_defocus[0] if len(cycle_defocus) >= 1 else np.nan
            cycle2 = cycle_defocus[1] if len(cycle_defocus) >= 2 else np.nan

            sem.L()
            ctf_defocus, ctf_res = run_ctffind()
            delta = af_defocus_used - ctf_defocus

            sem.Echo(
                f"target={target_defocus:.3f} um, corr={corr:.4f}, "
                f"AF(last of {autofocus_cycles})={af_defocus_used:.4f} um, "
                f"CTF={ctf_defocus:.4f} um, delta={delta:.4f} um, "
                f"drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
            )

            measurement_rows.append(
                {
                    "target_defocus_um": target_defocus,
                    "beam_tilt_correction": corr,
                    "autofocus_cycle_used": autofocus_cycles,
                    "autofocus_defocus_cycle1_um": cycle1,
                    "autofocus_defocus_cycle2_um": cycle2,
                    "autofocus_defocus_um": af_defocus_used,
                    "ctf_defocus_um": ctf_defocus,
                    "delta_um": delta,
                    "ctf_resolution_A": ctf_res,
                    "drift_speed_x_nm_per_s": speed_x,
                    "drift_speed_y_nm_per_s": speed_y,
                }
            )
            corr_vals.append(corr)
            delta_vals.append(delta)

        slope, intercept, best_correction = fit_best_correction(corr_vals, delta_vals)
        fitted_corrections.append(best_correction)
        summary_rows.append(
            {
                "target_defocus_um": target_defocus,
                "fitted_correction": best_correction,
                "fit_slope": slope,
                "fit_intercept": intercept,
            }
        )
        grouped_results.append(
            {
                "target_defocus": target_defocus,
                "corrections": corr_vals,
                "deltas": delta_vals,
                "slope": slope,
                "intercept": intercept,
                "best_correction": best_correction,
            }
        )
        sem.Echo(
            f"Fitted correction for defocus {target_defocus:.3f} um: {best_correction:.6f}"
        )

    finite_corr = [v for v in fitted_corrections if np.isfinite(v)]
    mean_correction = float(np.mean(finite_corr)) if finite_corr else float("nan")

    save_measurements_csv(csv_measurements, measurement_rows)
    save_summary_csv(csv_summary, summary_rows, mean_correction)
    make_plot(plot_file, grouped_results)

    sem.Echo("------------------------------------------------")
    sem.Echo("Fitted beam tilt correction values by target defocus:")
    for row in summary_rows:
        sem.Echo(
            f"  defocus {float(row['target_defocus_um']):.3f} um -> "
            f"correction {float(row['fitted_correction']):.6f}"
        )
    sem.Echo(f"Mean beam tilt correction: {mean_correction:.6f}")
    sem.Echo(f"Saved measurements CSV: {csv_measurements}")
    sem.Echo(f"Saved summary CSV: {csv_summary}")
    sem.Echo(f"Saved fit plot: {plot_file}")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
