#!Python
# ===================================================================
# ScriptName     LAFIS calibration
# Purpose:       calibrate laser alignment aberration free image shift
# Author:        Anchi Cheng, Misha Kopylov
# ===================================================================

import csv
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

############ SETTINGS ############
image_shift_delta_um = [(-5.0,0.0),]

xt_is_matrix = [[0.000312, -0.000352],[0.001120, 0.000281]]
sem.Echo("***Running X-Tilt LAFIS***")
sem.Echo("DIAGNOSTIC MESSAGE: Current image shift:")
sem.ReportImageShift()
sem.Echo("DIAGNOSTIC MESSAGE: Current beam tilt:")
sem.ReportBeamTilt()
sem.Echo("script starts here")

sem.SetImageShift(0,0)
#sem.SetXLensDeflector(2, -0.001909, -0.001152) #Reset to original X-tilt values
sem.Pause("Checkt X-tilt at no Image Shift")

def float_tuple(input_str_tuple):
   return (float(input_str_tuple[0]), float(input_str_tuple[1]))

def get_values():
   image_shift = sem.ReportImageShift()
   beam_tilt = sem.ReportBeamTilt()
   xtilt = sem.ReportXLensDeflector(2)
   return float_tuple(image_shift),float_tuple(beam_tilt), float_tuple(xtilt)

def calc_xt_is(xt0, is_delta):
   xt1 = [0.0,0.0]
   xt1[0] = xt0[0]+is_delta[0]*xt_is_matrix[0][0]+is_delta[1]*xt_is_matrix[1][0]
   xt1[1] = xt0[1]+is_delta[0]*xt_is_matrix[0][1]+is_delta[1]*xt_is_matrix[1][1]
   return xt1

is0, bt0, xt0 = get_values()
sem.R()

for delta in image_shift_delta_um:
   sem.SetImageShift(is0[0]+delta[0], is0[1]+delta[1])

   sem.AdjustBeamTiltforIS()
   is1, bt1, xt1 = get_values()
 
   sem.Echo("xt adjusted")
   xt = calc_xt_is(xt0, delta)
   sem.SetXLensDeflector(2, xt[0], xt[1])
   sem.R()
   sem.Pause("Adjust X-tilt")

   is2, bt2, xt2 = get_values()
   sem.R()
   for i in (0,1):
      if delta[i]:
         sem.Echo('-------%d-------' % i)
         sem.Echo('xt_delta %.6f, %.6f' % ((xt2[0]-xt0[0])/delta[i],(xt2[1]-xt0[1])/delta[i]))

sem.SetImageShift(is0[0],is0[1])
sem.SetXLensDeflector(2,xt0[0],xt0[1])
print(get_values())
