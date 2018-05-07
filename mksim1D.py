#!/usr/bin/python
# Creates a set of sensitivity experiments for the sea-ice ocean model
# F. Massonnet 2017

from sim1D import *
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

# Load namelist of experiments to perform
execfile("./sim1D_namelist.py")

n_exps = len(exps)

for j_exps in range(n_exps):
  print("Doing " + exps[j_exps][0])
  [state, t] = run_sim1D(k_ice        = exps[j_exps][1], \
                        alb_ice_mean = exps[j_exps][2], \
                        q_of_mean    = exps[j_exps][3], \
                        seed = 0,                       \
                        nt = 100 * 365)

  # Write NetCDF
  f = Dataset("./netcdfs/sim1D_" + exps[j_exps][0] + "_" + str(t[0].year) + str(t[0].month).zfill(2) + str(t[0].day).zfill(2) + "-" + str(t[-1].year) + str(t[-1].month).zfill(2) + str(t[-1].day).zfill(2)  + ".nc", mode = "w")

  f.Description = "Sea ice thickness from a 1-D sea ice - ocean model. Integration: " + exps[j_exps][0]
 
  time = f.createDimension("time", None)

  times = f.createVariable("time", np.int32, ("time",))
  times[:] = range(len(t))

  tmp = f.createVariable("sivol_Arctic_n80", np.float32, ("time",))
  # The volume is multiplied by the oceanic area north of 80N (http://www.bibmath.net/dico/index.php?action=affiche&quoi=./c/calotte.html)
  tmp[:] = state[0, :] * 2.0 * np.pi * 6371000.0 ** 2 * (1 - np.sin(80.0 / 360.0 * 2 * np.pi)) / 1e12
  tmp.units = "10^3 km3"
  del tmp

  tmp = f.createVariable("siarea_Arctic_n80", np.float32, ("time",))
  tmp[:] = state[1, :] * 2.0 * np.pi * 6371000.0 ** 2 * (1 - np.sin(80.0 / 360.0 * 2 * np.pi)) / 1e12         
  tmp.units = "10^6 km2"
  del tmp

  f.close()
