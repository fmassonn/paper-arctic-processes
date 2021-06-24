#!/usr/bin/python
#
# Script to make Fig. ED6 of the Feedback paper
# intended to highlight the link between process
# and variability. And plot PIOMAS
#
# Also produces statistics for Table 1

# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
exec(open("./config.py").read()) 


# Start the script
domain = "Arctic_n80" # Domain under consideration (north of X)
dirin = "./netcdfs/"   # where the pre-processed NetCDF data of volume and area are

# Start the script
# ----------------
# Create the figure
fig=plt.figure(figsize = (6, 10))

# PART 1 - Load and analyse CMIP5 data
# ------------------------------------
# Years defining the period to showcase (must be within the availability of data)
yearb = 1955
yeare = 2004
# Read CMIP5 namelist (in another file since this namelist is used by different scripts)
exec(open("./cmip5_namelist.py").read())

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
    nfb, _, _ = negative_seaice_feedback(volume, period = 12)
    # Melt
    pfb, _, _ = positive_seaice_feedback(volume, area, period = 12)
    # Amplitude seasonal cycle
    ampvol = np.mean([np.max(volume[t:t + 12]) - np.min(volume[t:t + 12]) for t in np.arange(0, len(volume), 12)])
    # Standard deviation March Volume
    stdv = np.std(volume[2::12])
    # Standard deviation September area
    stda = np.std(area[8::12])
    # Persistence
    vol_d = detrend(volume, order = 1, period = 12)
    pers = eftimescale(vol_d)[0]
    # March volume loss
    vchange = np.polyfit(range(len(volume[2::12])), volume[2::12], 1)[0] * (yeare - yearb)
    # September area loss
    achange = np.polyfit(range(len(area[8::12])), area[8::12], 1)[0] * (yeare - yearb)

    data[j_models].append([[nfb,    pfb ,   nfb,  pfb,     nfb,    pfb,  nfb,     pfb], \
                           [ampvol, ampvol, stdv, stda,    pers,   pers, vchange, achange]])

    fb_label = ["IFE [m/m]", "OWFE [m$^{-1}$]", \
                "IFE [m/m]", "OWFE [m$^{-1}$]", \
                "IFE [m/m]", "OWFE [m$^{-1}$]", \
                "IFE [m/m]", "OWFE [m$^{-1}$]", \
               ]

    va_label = ["Amplitude seasonal cycle\nice volume [10$^3$ km$^3$]", "Amplitude seasonal cycle\nice volume [10$^3$ km$^3$]", \
                "March ice-volume\n variability [10$^3$ km$^3$]",       "September ice-area\n variability [10$^6$ km$^2$]", \
                "Ice-volume\npersistence [months]",                      "Ice-volume\npersistence [months]", \
                "March ice-volume\nloss [10$^3$ km$^3$]",               "September ice-area\nloss [10$^6$ km$^2$]", \
               ]

    xlims = [
             [-1.0, 0.2], [-0.2, 1.0], \
             [-1.0, 0.2], [-0.2, 1.0], \
             [-1.0, 0.2], [-0.2, 1.0], \
             [-1.0, 0.2], [-0.2, 1.0], \
            ]

    ylims = [
             [-0.5, 10.0], [-0.5, 10.0], \
             [-0.1, 2.0], [-0.03, 0.6] , \
             [-5.0, 100.0], [-5.0, 100.0], \
             [-5.0, 0.25], [-1.5, 0.075] , \
            ]

    for jsubplot in range(8):  # For each pair of diagnostic
      plt.subplot(4, 2, jsubplot + 1)
      plt.scatter(data[j_models][j_members][0][jsubplot], data[j_models][j_members][1][jsubplot], 1, color = colors[j_models])
      plt.xlabel(fb_label[jsubplot])
      plt.xlim(xlims[jsubplot])
      plt.ylim(ylims[jsubplot])
      plt.ylabel(va_label[jsubplot])
      plt.plot((-1e9, 1e9), (0, 0), 'k', lw = 2)
      plt.plot((0, 0), (-1e9, 1e9), 'k', lw = 2)

      # If this is the last member, plot the ensemble mean diagnostics
      if j_members == n_members - 1:
        x = np.mean([data[j_models][jj][0][jsubplot] for jj in range(n_members)])
        y = np.mean([data[j_models][jj][1][jsubplot] for jj in range(n_members)])
        plt.scatter(x, y, 20 * n_members, color = colors[j_models], alpha = 0.5)
        plt.text(   x, y, str(j_models + 1).zfill(2), color = [0.5 * c for c in colors[j_models]], fontsize = 4, ha = 'center', va = 'center')
  
        # If this is the last model, return the statistics
        if j_models == n_models - 1:
          x = [np.mean([data[j_models][jj][0][jsubplot] for jj in range(n_members)]) for j_models in range(n_models)]
          y = [np.mean([data[j_models][jj][1][jsubplot] for jj in range(n_members)]) for j_models in range(n_models)]
          r = np.corrcoef(x, y)[0, 1]
          tstat = r / np.sqrt((1 - r ** 2) / (n_models - 2))  # The t-statistic.
                                                              # Under the null hypothesis of no correlation,
                                                              # tstat follows a student's law with  N - 2 dof.
          pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), n_models - 2)

          plt.title("R = " + str(np.round(r, 2)) + ", P = " + '{0:1.2e}'.format(pval))
      if j_models == 0 and j_members == 0:
        plt.grid()
        jsubplot += 1

# Get the PIOMAS estimates
# Add PIOMAS estimates
yearb = 1979
yeare = 2015
ncload(dirin + "piomas_" + str(yearb) + "01-" + str(yeare) + "12.nc", ["sivol_" + domain, "siarea_" + domain])
volume = eval("sivol_" + domain)
area = eval("siarea_" + domain)

# Amplitude seasonal cycle
ampvol = np.mean([np.max(volume[t:t + 12]) - np.min(volume[t:t + 12]) for t in np.arange(0, len(volume), 12)])
# Standard deviation March Volume
stdv = np.std(volume[2::12])
# Standard deviation September area
stda = np.std(area[8::12])
# Persistence
vol_d = detrend(volume, order = 1, period = 12)
pers = eftimescale(vol_d)[0]
# March volume loss
vchange = np.polyfit(range(len(volume[2::12])), volume[2::12], 1)[0] * (yeare - yearb)
# September area loss
achange = np.polyfit(range(len(area[8::12])), area[8::12], 1)[0] * (yeare - yearb)

# Evaluate the processes
nf, stats_nf, _ = negative_seaice_feedback(volume, period = 12)
pf, stats_pf, _ = positive_seaice_feedback(volume, area, period = 12)

plt.subplot(4, 2, 1); plt.scatter(nf, ampvol, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 2); plt.scatter(pf, ampvol, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 3); plt.scatter(nf, stdv, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 4); plt.scatter(pf, stda, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 5); plt.scatter(nf, pers, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 6); plt.scatter(pf, pers, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 7); plt.scatter(nf, vchange, 300, marker = '*', color = 'black')
plt.subplot(4, 2, 8); plt.scatter(pf, achange, 300, marker = '*', color = 'black')
       
plt.tight_layout()
plt.savefig("./figED6.png", dpi = 500)
