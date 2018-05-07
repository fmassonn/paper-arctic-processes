#!/usr/bin/python
#
# Script to make ED Fig. 11 of the process paper
# intended to highlight the (non-) sensitivity of
# results to detrending 
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
# Create the figure
fig=plt.figure(figsize = (4, 6))
gs=GridSpec(2, 1) # 2 rows, 2 columns

ax = []
ax.append(fig.add_subplot(gs[0, 0]))
ax.append(fig.add_subplot(gs[1, 0]))

# Load and analyse CMIP5 data
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
  # 1. Get information about the model
  model = info[j_models][0]
  yb    = info[j_models][5]
  ye    = info[j_models][6]

  members = info[j_models][7]
  n_members = len(members)

  print("Loading " + model + "(" + str(n_members) + " members)")

  data.append(list()) # One per member
  for j_members in range(n_members):
    # 2. Read volume, area and temperature
    ncload(dirin + "/seaice_" + model + "_historical_r" + \
           str(members[j_members]) + "i1p1_" + str(yb) +  \
           "01-" + str(ye) + "12.nc", ["sivol_" + domain, \
           "siarea_" + domain])

    volume = eval("sivol_" + domain)
    # Subset to the years of interest
    volume = volume[(yearb - yb) * 12:(yeare - yb) * 12 + 12]
  
    area = eval("siarea_" + domain)
    # Subset to the years of interest
    area = area[(yearb - yb) * 12:(yeare - yb) * 12 + 12]


    # Mean state
    vmean = np.mean(volume)

    # Growth
    nfb, _, _ = negative_seaice_feedback(volume, period = 12, order = 0)
    # Melt
    pfb, _, _ = positive_seaice_feedback(volume, area, period = 12, order = 0)

    data[j_models].append([vmean, nfb, pfb])


# Scatter plots
# -------------
comp = info[j_models][9] - 1 # Complexity index
fit = list()
std = list()

for j_models in range(n_models):
  comp = info[j_models][9] - 1 # Complexity index retrieved from namelist
  members = info[j_models][7]
  n_members = len(members)

  # Plot mean of diagnostics
  x  = np.mean([i[0] for i in data[j_models]])
  y1 = np.mean([i[1] for i in data[j_models]])
  y2 = np.mean([i[2] for i in data[j_models]])
 
  ax[0].scatter(x, y1, 20 * n_members, color = colors[j_models], alpha = 0.5)
  ax[0].text(   x, y1, str(j_models + 1).zfill(2), color = [c / 2.0 for c in colors[j_models]], fontsize = 4, ha = 'center', va = 'center')
  ax[1].scatter(x, y2, 20 * n_members, color = colors[j_models], alpha = 0.5)
  ax[1].text(   x, y2, str(j_models + 1).zfill(2), color = [c / 2.0 for c in colors[j_models]], fontsize = 4, ha = 'center', va = 'center')
  del x, y1, y2

  # Plot diagnostics per member
  for j_members in range(n_members):
    x = data[j_models][j_members][0]
    y1 = data[j_models][j_members][1]
    y2 = data[j_models][j_members][2]
    ax[0].scatter(x, y1, 0.5, color = colors[j_models])
    ax[1].scatter(x, y2, 0.5, color = colors[j_models])


#ax[0].set_xlabel("Annual-mean volume [10$^3$ km$^3$]")
ax[0].set_ylabel("IFE [m/m]")
ax[0].plot((-1e9, 1e9), (0, 0), 'k-', lw = 3)
ax[0].plot((0, 0), (-1e9, 1e9), 'k-', lw = 3)
ax[0].set_xlim(-1.0, 25.0)
ax[0].set_ylim(-1.0, 0.2 )
ax[0].set_title("44 CMIP5 models (146 members)")
ax[0].text(0.7, 0.1, "a", fontweight = "bold", fontsize = 16)
ax[0].grid()

ax[1].set_xlabel("Annual mean sea-ice volume [10$^3$ km$^3$]")
ax[1].set_ylabel("OWFE [m$^{-1}$]")
ax[1].plot((-1e9, 1e9), (0, 0), 'k-', lw = 3)
ax[1].plot((0, 0), (-1e9, 1e9), 'k-', lw = 3)
ax[1].set_xlim(-1.0, 25.0)
ax[1].set_ylim(-0.2, 1.0 )
ax[1].text(0.7, 0.9, "b", fontweight = "bold", fontsize = 16)
ax[1].grid()
  

plt.tight_layout()
plt.savefig("./figED11.png", dpi = 500)
