#!/usr/bin/python
#
# Script to make ED Figs. 8, 9, 12 and 13 of the process paper
# intended to show basic time series of volume, area, and the
# scatter plots leading to the process
# Loads the first member only
#

# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


# Start the script
domain = "Arctic_n80" # Domain under consideration (north of X)
dirin = "./netcdfs/"   # where the pre-processed NetCDF data of volume and area are

# Start the script
# ----------------
# Create the figures
[plt.figure("fig" + str(j + 1), figsize = (12, 15)) for j in range(4)]

# PART 1 - Load and analyse CMIP5 data
# ------------------------------------
# Years defining the period to showcase (must be within the availability of data)
yearb = 1955
yeare = 2004
# Read CMIP5 namelist (in another file since this namelist is used by different scripts)
execfile("./cmip5_namelist.py")

n_models = len(info)
n_years = yeare - yearb + 1

# Plotting colors
np.random.seed(1) # to fix the colors of CMIP5
colors = [np.random.random(3) for j in range(n_models)]

# Load the data and compute the diagnostics
data = list() # Will contain as many items as there are models
              # Each item will contain as many sub-items as there are members
nt = (yeare - yearb) * 12 + 1 # Number of time steps (months) in the analysis

for j_models in range(n_models):
  for j in range(4):
    plt.figure("fig" + str(j + 1))
    plt.subplot(8, 6, j_models + 1)
  
  # 1. Get information about the model
  model = info[j_models][0]
  yb    = info[j_models][5]
  ye    = info[j_models][6]

  members = info[j_models][7]
  n_members = len(members)

  print("Loading " + model)

  data.append(list()) # One per member
  for j_members in [0]:
    # 2. Read volume, area and temperature
    ncload(dirin + "/seaice_" + model + "_historical_r" + \
           str(members[j_members]) + "i1p1_" + str(yb) +  \
           "01-" + str(ye) + "12.nc", ["sivol_" + domain, \
           "siarea_" + domain])

    volume = eval("sivol_" + domain)
    # Subset to the years of interest
    volume_subset = volume[(yearb - yb) * 12:(yeare - yb) * 12 + 12]
  
    area = eval("siarea_" + domain)
    # Subset to the years of interest
    area_subset = area[(yearb - yb) * 12:(yeare - yb) * 12 + 12]
    

    plt.figure("fig1")
    plt.plot(np.arange(yb, ye + 1, 1.0 / 12), volume, color = colors[j_models])
    plt.grid()
    plt.ylim(0, 25)
    plt.xticks(np.arange(1850, 2005, 50))
    plt.title(model)

    plt.figure("fig2")
    plt.plot(np.arange(yb, ye + 1, 1.0 / 12), area, color = colors[j_models])
    plt.grid()
    plt.ylim(0, 5)
    plt.xticks(np.arange(1850, 2005, 50))
    plt.title(model)

    # Growth
    nf, stats, series = negative_seaice_feedback(volume_subset, period = 12)

    plt.figure("fig3")
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.plot((-1e9, 1e9), (0, 0), 'k')
    plt.plot((0, 0), (-1e9, 1e9), 'k')
    plt.plot((-10, 10), (-10 * nf, 10 * nf), color = colors[j_models])
    plt.grid()
    plt.scatter(series[0], series[1], 8, color = colors[j_models], marker = ',')
    plt.title(model + "\ny = " + str(np.round(nf, 2)) + " x")

    # Melt
    pf, stats, series = positive_seaice_feedback(volume_subset, area_subset, period = 12)

    plt.figure("fig4")
    plt.xlim(-2.5, 2.5)
    plt.ylim(-1.0, 1.0)
    plt.plot((-1e9, 1e9), (0, 0), 'k')
    plt.plot((0, 0), (-1e9, 1e9), 'k')
    plt.plot((-10, 10), (-10 * pf, 10 * pf), color = colors[j_models])
    plt.grid()
    plt.scatter(series[0], series[1], 8, color = colors[j_models], marker = ',')
    plt.title(model + "\ny = " + str(np.round(pf, 2)) + " x")



figindex = [8, 9, 12, 13]
for j in range(4):
  plt.figure("fig" + str(j + 1))
  plt.tight_layout()
  plt.savefig("./figED" + str(figindex[j]) + ".png", dpi = 300)
  plt.close("fig" + str(j + 1))

plt.clf()



