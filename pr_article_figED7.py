#!/usr/bin/python
#
# Script to make ED Fig. 7 of the process paper
# Toy model in reference case

# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
exec(open("./config.py").read())

from sim1D import *

[state, t] = run_sim1D(seed = 1)

offset = [0.0, 0.0, -273.15, -273.15]
scalef = [1.0, 100.0, 1.0, 1.0]
lab = ["Volume [m]", "Concentration [%]", "Ice temperature [$^\circ$C]", "Mixed layer\ntemperature [$^\circ$C]"]
plt.figure(figsize = (6, 10))
for j in range(4):
  plt.subplot(4, 1, j + 1)
  plt.plot(t, scalef[j] * state[j, :] + offset[j], lw = 2)
  plt.xlim(t[0], t[-1])
  plt.grid()
  plt.ylabel(lab[j])

plt.tight_layout()
plt.savefig("./figED7.png", dpi = 300)
plt.clf()


