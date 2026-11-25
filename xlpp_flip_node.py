#!Python
import ronchi_cal
import numpy as np

ronchi_cal.ronchi_sem_lib.checkRonchigramSetup()
corr_scale = 1e-5
sem=ronchi_cal.sem
m=ronchi_cal.ronchi_sem_lib.ronchiCorrMatrix

def flip(sign=1):
    x0=sem.ReportXLensDeflector(2)
    corr = np.asarray(m, dtype=float)*corr_scale
    xt = x0 + np.pi * sign * corr[0] + np.pi * sign * corr[1]
    sem.SetXLensDeflector(2,xt[0],xt[1])

import sys
if len(sys.argv) != 2:
    print("usage: python xlpp_flip_node.py <scale to multiply calibrated pi radian phase shift to add (1 or -1)>")
    sys.exit(1)
r = int(sys.argv[1])

if r in (-1,1):
    print(f'shift xt by ({r:d}*pi,{r:d}*pi)')
    flip(r)
