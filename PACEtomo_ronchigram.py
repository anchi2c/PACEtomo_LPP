#!Python
# Ronchigram FFT analysis and Xtilt correction for PACEtomo

import numpy as np


def bin_image(image, binning=32):
    bins = [image[i:(image.shape[0] // binning * binning):binning, j:(image.shape[1] // binning * binning):binning]
            for i in range(binning) for j in range(binning)]
    return np.sum(bins, axis=0)


def find_fourier_centered(image, padded_size=4096):
    fourier = np.fft.fftshift(np.fft.fft2(image - np.mean(image), s=(padded_size, padded_size)))
    center = np.array(image.shape) / 2
    grid_y, grid_x = np.indices((padded_size, padded_size), dtype=float)
    shifted_y = grid_y - padded_size / 2
    shifted_x = grid_x - padded_size / 2
    correction_phase_x = np.exp(2j * np.pi * shifted_x / padded_size * center[1])
    correction_phase_y = np.exp(2j * np.pi * shifted_y / padded_size * center[0])
    correction_phase = correction_phase_x * correction_phase_y
    return fourier * correction_phase


def report_angles(ks, start_angle=-135):
    return np.mod(np.arctan2(-ks[:, 0], ks[:, 1]) * 180 / np.pi - start_angle, 360) + start_angle


def find_peaks(fourier, radius=100, npeaks=4):
    fourier_abs = np.abs(fourier).copy()
    peak_locations = []
    phases = []
    size = np.shape(fourier_abs)[0]
    peak_coords = [size // 2, size // 2]
    fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0

    for _ in range(npeaks):
        max_idx = np.argmax(fourier_abs)
        peak_coords = np.unravel_index(max_idx, fourier_abs.shape)
        peak_locations.append(peak_coords)
        fourier_abs[(peak_coords[0] - radius):(peak_coords[0] + radius),
                    (peak_coords[1] - radius):(peak_coords[1] + radius)] = 0
        phases.append(np.angle(fourier[peak_coords]))

    return np.array(peak_locations) - size / 2, np.array(phases)


def find_ks_phases(corrected_fourier, pixel_size_um, npeaks=2, radius=100, binning=1, fourier_size=None):
    if fourier_size is None:
        fourier_size = corrected_fourier.shape[0]
    peaks, phases = find_peaks(corrected_fourier, radius=radius, npeaks=npeaks * 2)
    ordering = np.argsort(report_angles(peaks, start_angle=-135))
    peaks = peaks[ordering][:npeaks]
    phases = phases[ordering][:npeaks]
    ks = peaks / fourier_size * 1 / (pixel_size_um * binning)
    return ks, phases


def analyze_ronchigram(image, pixel_size_um, binning, target_phase_a, target_phase_b,
                       correct_ks, peak_radius=100, corr_matrix=None, corr_scale=1e-5):
    """Analyze ronchigram image and compute XLensDeflector correction deltas."""
    if corr_matrix is None:
        corr_matrix = [[0.212, 1.28], [1.22, -0.243]]

    correct_ks = np.asarray(correct_ks, dtype=float)
    binned = bin_image(np.asarray(image), binning=binning)
    image_fft = find_fourier_centered(binned)
    ks, phases = find_ks_phases(image_fft, pixel_size_um * binning, npeaks=2, radius=peak_radius, binning=1,
                                fourier_size=image_fft.shape[0])

    phase_err_a = np.mod(phases[0] - target_phase_a + np.pi, 2 * np.pi) - np.pi
    phase_err_b = np.mod(phases[1] - target_phase_b + np.pi, 2 * np.pi) - np.pi

    corr = np.asarray(corr_matrix, dtype=float) * corr_scale
    correction_x = phase_err_a * corr[0, 0] + phase_err_b * corr[0, 1]
    correction_y = phase_err_a * corr[1, 0] + phase_err_b * corr[1, 1]

    return {
        "ks": ks,
        "phases": phases,
        "ks_error": ks - correct_ks,
        "phase_err_a": phase_err_a,
        "phase_err_b": phase_err_b,
        "correction_x": correction_x,
        "correction_y": correction_y,
    }


def apply_laser_correction(sem_module, correction_x, correction_y, lens_index=2):
    xtX, xtY = sem_module.ReportXLensDeflector(lens_index)
    new_x = xtX + correction_x
    new_y = xtY + correction_y
    sem_module.SetXLensDeflector(lens_index, new_x, new_y)
    return new_x, new_y
