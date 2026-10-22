#!Python
import ronchi_cal

ronchi_cal.ronchi_sem_lib.checkRonchigramSetup()
pixel_size_angs = ronchi_cal.sem.ReportCurrentPixelSize('T')
xt_tilt = 1 #(scaled at 1e-5 rad)
corr_scale = 1e-5
image_binning = 2
ronchi_binning = 32
# The same value has always been used.
ronchi_fixed_pixel_size_um = 0.98e-4 * 2

if __name__=='__main__':
    # ronchiStartC3 calibration
    new_ronchi_start_c3 = ronchi_cal.calibrate_ronchigram_start_c3(ronchi_fixed_pixel_size_um, ronchi_binning,
                       xt_tilt, peak_radius=100, corr_scale=corr_scale,
                       c3_correction_factor=20 / 6.85)
    print('new StartC3Offset',new_ronchi_start_c3)
