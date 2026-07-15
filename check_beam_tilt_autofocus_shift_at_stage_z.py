#!Python
# ===================================================================
# ScriptName     test_beam_tilt
# Purpose:       plot beam tilt image shift at different stage z height
#                compared with CtfFind defocus.
# Author:        Generated for PACEtomo workflow
# ===================================================================

import csv
import time
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

############ SETTINGS ############

# Beam tilt used for autofocus measurement [mrad]
sc=1.7335
tilt_angle_mrad_values = [-10.0,10.0]
tilt_angle_mrad_values = list(map((lambda x: x*sc), tilt_angle_mrad_values))

# Correction factors to beam tilt scale in conversion to defocus
correction = 1

# Target stage z values [microns] to calibrate at
# 1 um backlash applied before the first
target_stage_z_values = [1.5,1.0,0.5,0.0,-0.5,-1,-1.5]

# CTF search range [microns]
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

# Defocus convergence settings for setting each target defocus
target_stage_z_tolerance_um = 0.01
max_defocus_adjust_iterations = 5

# Number of autofocus cycles per test.
# Delta is computed using the last cycle's autofocus value.
autofocus_cycles = 1

# Output file names (written in current SerialEM working directory)
csv_measurements = "beam_tilt_autofocus_measurements.csv"
csv_summary = "beam_tilt_autofocus_summary.csv"
plot_file = "beam_tilt_autofocus_fits.png"

########## END SETTINGS ##########

xyz0 = sem.ReportStageXYZ()

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


def set_target_stage_z(target_stage_z):
    """Adjust focus to target defocus and verify with CtfFind."""
    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    sem.MoveStageTo(xyz0[0],xyz0[1],xyz0[2]+target_stage_z)
    time.sleep(5)
    sem.ReportStageXYZ()
    current_defocus = np.nan
    for attempt in range(1):
        sem.L()
        current_defocus, _ = run_ctffind()
        error = target_stage_z - current_defocus
        sem.Echo(
            f"Target defocus check {attempt}/{max_defocus_adjust_iterations}: "
            f"current={current_defocus:.3f} um, target={target_stage_z:.3f} um, "
            f"error={error:.3f} um"
        )
        if True:
            return current_defocus
    sem.Echo(
        "WARNING: Target defocus not reached within tolerance after max iterations."
    )
    return current_defocus


def run_autofocus_trial(tilt_angle_mrad):
    """
    Execute one autofocus-like measurement with a given beam tilt_angle.
    Returns:
      defocus_measured [microns], speed_x [nm/s], speed_y [nm/s]
    """
    beam_tilt = sem.ReportBeamTilt()
    tilt_x_orig = float(beam_tilt[0])
    tilt_y_orig = float(beam_tilt[1])
    tilt_x_plus = tilt_x_orig + correction * tilt_angle_mrad

    pixel_size_binned = float(sem.ReportCurrentPixelSize("R"))
    binning = float(sem.ReportBinning("R"))
    pixel_size_unbinned = pixel_size_binned / binning

    # Positive tilt
    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    bt = sem.ReportBeamTilt()
    sem.Echo(f"ActuralBeamTilt: {float(bt[0]):.5f} ")
    sem.F() #acquire image with Focus preset
    sem.ResetClock()
    sem.Copy("A", "L")

    # Zero tilt
    sem.SetBeamTilt(tilt_x_orig, tilt_y_orig)
    sem.F()
    sem.AlignTo("L", 1)
    align_shift_1 = sem.ReportAlignShift()
    disp_x1_px = float(align_shift_1[0])
    disp_y1_px = float(align_shift_1[1])

    # Positive tilt again to correct for drift
    sem.SetBeamTilt(tilt_x_plus, tilt_y_orig)
    sem.F()
    sem.AlignTo("L", 1)
    align_shift_2 = sem.ReportAlignShift()
    disp_x2_px = float(align_shift_2[0])
    disp_y2_px = float(align_shift_2[1])
    elapsed = float(sem.ReportClock())

    # Align to origin
    sem.SetBeamTilt(tilt_x_orig, tilt_y_orig)

    drift_x = disp_x2_px * pixel_size_unbinned
    drift_y = disp_y2_px * pixel_size_unbinned
    speed_x = drift_x / elapsed if elapsed > 0 else 0.0
    speed_y = drift_y / elapsed if elapsed > 0 else 0.0

    displacement_from_tilt_x = (disp_x1_px - disp_x2_px / 2.0) * pixel_size_unbinned
    displacement_from_tilt_y = (disp_y1_px - disp_y2_px / 2.0) * pixel_size_unbinned
    displacement = np.sqrt(
        displacement_from_tilt_x * displacement_from_tilt_x
        + displacement_from_tilt_y * displacement_from_tilt_y
    )

    if displacement_from_tilt_x == 0:
        sign = 1.0
    else:
        sign = displacement_from_tilt_x / abs(displacement_from_tilt_x)

    defocus_measured = -1.0 * sign * displacement

    return defocus_measured, speed_x, speed_y

def fit_auto_focus(corrections, af_vals):
    """Linear fit of delta vs correction; returns slope, intercept, root."""
    result, residuals, rank, singular_values, rcond = np.polyfit(corrections, af_vals, 1, full=True)
    slope = result[0]
    intercept = result[1]
    return slope, intercept, residuals

def save_measurements_csv(path, rows):
    fields = [
        "target_stage_z_um",
        "tilt_angle_mrad",
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


def save_summary_csv(path, rows):
    fields = [
        "target_stage_z_um",
        "fit_slope",
        "fit_intercept",
        "fit_residuals",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow(
            {
                "target_stage_z_um": "MEAN",
                "fit_slope": "",
                "fit_intercept": "",
                "fit_residuals": "",
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
        print(res)
        x = np.array(res["tilt_angle_mrad_values"], dtype=float)
        y = np.array(res["af_values"], dtype=float)
        slope = res["slope"]
        intercept = res["intercept"]
        best = res["residuals"]

        ax.scatter(x, y, label="measurements")
        xfit = np.linspace(np.min(x), np.max(x), 200)
        ax.plot(xfit, slope * xfit + intercept, label="linear fit")
        ax.axhline(0, color="gray", linestyle="--", linewidth=1)
        if np.isfinite(slope):
            ax.axvline(slope, color="red", linestyle=":", linewidth=1)
            title = f"Stage Z {res['target_stage_z']:.2f} um, slope={slope:.4f}, intercept={intercept:.4f}"
        else:
            title = f"Stage Z {res['target_stage_z']:.2f} um, slope=NaN"
        ax.set_title(title)
        ax.set_xlabel("Beam tilt (mrad)")
        ax.set_ylabel("Autofocus value (?)")
        ax.legend()

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    sem.Echo("##### Beam tilt calibration test #####")
    sem.Echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    sem.Echo(f"Testing beam tilts: {list(np.array(tilt_angle_mrad_values, dtype=float))}")
    sem.Echo(f"Target stage values from eucentric: {target_stage_z_values}")

    measurement_rows = []
    summary_rows = []
    grouped_results = []

    set_target_stage_z(target_stage_z_values[0]+1.0)
 
    for target_stage_z in target_stage_z_values:
        sem.Echo("------------------------------------------------")
        sem.Echo(f"Measuring target defocus for stage z {target_stage_z:.3f} um")
        reached_defocus = set_target_stage_z(target_stage_z)
        sem.Echo(f"CtffindDefocus: {reached_defocus:.3f} um")

        tilt_vals = []
        af_vals = []
        for tilt in tilt_angle_mrad_values:
            tilt = float(tilt)
            cycle_defocus = []
            speed_x = 0.0
            speed_y = 0.0
            for _ in range(autofocus_cycles):
                af_defocus, speed_x, speed_y = run_autofocus_trial(tilt)
                cycle_defocus.append(af_defocus)
            af_defocus_used = cycle_defocus[-1]
            cycle1 = cycle_defocus[0] if len(cycle_defocus) >= 1 else np.nan
            cycle2 = cycle_defocus[1] if len(cycle_defocus) >= 2 else np.nan

            sem.L()
            ctf_defocus, ctf_res = run_ctffind()
            delta = af_defocus_used - ctf_defocus

            sem.Echo(
                f"target={target_stage_z:.3f} um, tilt={tilt:.4f}, "
                f"AF(last of {autofocus_cycles})={af_defocus_used:.4f} um, "
                f"CTF={ctf_defocus:.4f} um, delta={delta:.4f} um, "
                f"drift=({speed_x:.3f}, {speed_y:.3f}) nm/s"
            )

            measurement_rows.append(
                {
                    "target_stage_z_um": target_stage_z,
                    "tilt_angle_mrad": tilt,
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
            tilt_vals.append(tilt)
            af_vals.append(af_defocus_used)

        slope, intercept, residuals = fit_auto_focus(tilt_vals, af_vals)
        summary_rows.append(
            {
                "target_stage_z_um": target_stage_z,
                "fit_slope": slope,
                "fit_intercept": intercept,
                "fit_residuals": residuals,
            }
        )
        grouped_results.append(
            {
                "target_stage_z": target_stage_z,
                "tilt_angle_mrad_values": tilt_vals,
                "af_values": af_vals,
                "slope": slope,
                "intercept": intercept,
                "residuals": residuals,
            }
        )
        sem.Echo(
            f"Fitted slope for defocus {target_stage_z:.3f} um: {slope:.6f}, "
            f"Fitted intercept for defocus {target_stage_z:.3f} um: {intercept:.6f}"
        )

    sem.MoveStageTo(xyz0[0],xyz0[1],xyz0[2])

    save_measurements_csv(csv_measurements, measurement_rows)
    save_summary_csv(csv_summary, summary_rows)
    make_plot(plot_file, grouped_results)

    sem.Echo("------------------------------------------------")
    sem.Echo("Fitted beam tilt correction values by target defocus:")
    for row in summary_rows:
        sem.Echo(
            f"  defocus {float(row['target_stage_z_um']):.3f} um -> "
        )
    sem.Echo(f"Saved measurements CSV: {csv_measurements}")
    sem.Echo(f"Saved summary CSV: {csv_summary}")
    sem.Echo(f"Saved fit plot: {plot_file}")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
