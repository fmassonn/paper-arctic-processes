#!/usr/bin/python
#
# Script to make Fig. 4 of the Process paper
# intended to explain the implications
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


domain = "Arctic" # Domain under consideration (north of X)
dirin = "./netcdfs/"   # where the pre-processed NetCDF data of volume and area are


# Create the figure
fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (10, 7))

fig.subplots_adjust(wspace=.5)

# Load the CESM-LENS data (1920-2100)
# -----------------------------------
n_members = 35

volume = np.empty(((2100 - 1920 + 1) * 12, n_members))
area   = np.empty(((2100 - 1920 + 1) * 12, n_members))
nt = (2100 - 1920 + 1) * 12 # Length of vector
volume[:] = np.nan
area[:] = np.nan

for j_members in np.arange(1, n_members + 1):
  print(j_members)
  if j_members == 1:
    ncload(dirin + "b.e11.B20TRC5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".185001-200512.nc", ["sivol_" + domain, "siarea_" + domain])
    volume1 = eval("sivol_" + domain + "[(1920 - 1850) * 12:]") # subsetting to adjust to other members, that start in 1920
    area1   = eval("siarea_" + domain + "[(1920 - 1850) * 12:]") 
    ncload(dirin + "b.e11.BRCP85C5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".200601-210012.nc", ["sivol_" + domain, "siarea_" + domain])
    volume2 = eval("sivol_" + domain)
    area2   = eval("siarea_" + domain)

  else:
    ncload(dirin + "b.e11.B20TRC5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".192001-200512.nc", ["sivol_" + domain, "siarea_" + domain])
    volume1 = eval("sivol_" + domain)
    area1   = eval("siarea_" + domain)
    ncload(dirin + "b.e11.BRCP85C5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".200601-210012.nc", ["sivol_" + domain, "siarea_" + domain])
    volume2 = eval("sivol_" + domain)
    area2   = eval("siarea_" + domain)

  volume[:, j_members - 1] = np.concatenate((volume1, volume2), 0)
  area[:, j_members - 1]   = np.concatenate((area1  , area2)  , 0)

# 1. Amplitude of seasonal cycle as function of mean state
# --------------------------------------------------------
# Hack to create colorbar
Z = [[0,0],[0,0]]
levels = range(1920 + 10, 2080 + 10 + 1, 1)
CS3 = ax[0][0].contourf(Z, levels, cmap = plt.cm.YlOrBr)
diag = (list(), list(), list(), list(), list()) # Will contain the mean (1) and variability (4)
                                                    # diagnostics for each member and each period
for j_members in range(n_members):
  print(j_members)
  [diag[i].append(list()) for i in range(5)] # Create new list to store all time-interval diags

  for year in np.arange(1920, 2080, 1):
    t1 = (year - 1920) * 12
    t2 = (year + 20 - 1920) * 12 + 12

    # Mean state
    vmean = np.mean(volume[t1:t2, j_members])
    diag[0][j_members].append(vmean)

    # Annual-mean volume
    vanmean = np.array([np.mean(volume[t:t + 12, j_members]) for t in np.arange(t1, t2, 12)])

    # Mean amplitude of seasonal cycle
    mamp = np.mean([np.max(volume[t:t + 12, j_members]) - np.min(volume[t:t + 12, j_members]) for t in np.arange(t1, t2, 12) ])
    diag[1][j_members].append(mamp)

    # Variability
    std  = np.std(vanmean)
    diag[2][j_members].append(std)

    # Persistence
    vol_d = detrend(volume[t1:t2, j_members], period = 12, order = 1)
    pers = eftimescale(vol_d)[0]
    diag[3][j_members].append(pers)

    # Trends
    tre  = np.polyfit(range(len(vanmean)), vanmean, 1)[0]
    diag[4][j_members].append(tre)

    # Record the diagnostic for averaging later
    ax[0][0].scatter(vmean, mamp, 1, color = plt.cm.YlOrBr(int(1.0 * (year + 10 - 1920) / 181 * 255)))
    ax[0][1].scatter(vmean, std , 1, color = plt.cm.YlOrBr(int(1.0 * (year + 10 - 1920) / 181 * 255)))
    ax[1][0].scatter(vmean, pers, 1, color = plt.cm.YlOrBr(int(1.0 * (year + 10 - 1920) / 181 * 255)))
    ax[1][1].scatter(vmean, tre , 1, color = plt.cm.YlOrBr(int(1.0 * (year + 10 - 1920) / 181 * 255)))

    ax[0][0].set_ylim(0, 20.0)
    ax[0][1].set_ylim(0, 5)
    ax[1][0].set_ylim(0, 40)
    ax[1][1].set_ylim(-1.0, 0.5)
    for a in ax.flat:
      a.plot((-1e9, 1e9), (0, 0), 'k--')
    for a in ax.flat:
      a.set_xlim(0, 44)
   
    del vmean, vanmean, mamp, std, vol_d, pers, tre

# Add reanalysis estimates
del volume, area

ncload(dirin + "piomas_197901-201512.nc", ["sivol_" + domain])
volume = eval("sivol_" + domain)
# 1. Compute the diagnostics on the full (1979-2015) period
vmean = np.mean(volume)
mamp = np.mean([np.max(volume[t:t + 12]) - np.min(volume[t:t + 12]) for t in np.arange(0, len(volume), 12)])
ax[0][0].scatter(vmean, mamp, 40, marker = 'x', color = [0.0, 0.0, 0.0])
vanmean = [np.mean(volume[t:t + 12]) for t in np.arange(0, len(volume), 12)]
vanstd = np.std(vanmean)
ax[0][1].scatter(vmean, vanstd, 40, marker = 'x', color = [0.0, 0.0, 0.0])
vol_d = detrend(volume, period = 12, order = 1)
pers = eftimescale(vol_d)[0]
ax[1][0].scatter(vmean, pers, 40, marker = 'x', color = [0.0, 0.0, 0.0])
tre  = np.polyfit(range(len(vanmean)), vanmean, 1)[0]
ax[1][1].scatter(vmean, tre, 40, marker = 'x', color = [0.0, 0.0, 0.0])

# 2. Compute the diagnostics on sub-periods of 20-yr as a measure
#    of uncertainty
vmean_r = list()
mamp_r  = list()
vanstd_r= list()
pers_r  = list()
tre_r   = list()

year = 1979
while year + 20 <= 2015:
  t1 = (year - 1979) * 12
  t2 = (year + 20 - 1979) * 12 + 12
  
  vmean_r.append(np.mean(volume[t1:t2]))
  vanmean = np.array([np.mean(volume[t:t + 12]) for t in np.arange(t1, t2, 12)])

  mamp_r.append(np.mean([np.max(volume[t:t + 12]) - np.min(volume[t:t + 12]) for t in np.arange(t1, t2, 12)]))

  vanstd_r.append(np.std(vanmean))

  vol_d = detrend(volume[t1:t2], period = 12, order = 1)
  pers_r.append(eftimescale(vol_d)[0])

  tre_r.append(np.polyfit(range(len(vanmean)), vanmean, 1)[0]) 

  year += 1

ax[0][0].errorbar(vmean, mamp,   xerr = np.std(vmean_r), yerr = np.std(mamp_r),   color = [0.0, 0.0, 0.0])
ax[0][0].annotate('Reanalysis', xy=(vmean - 0.5, mamp - 0.3), xytext=(25, 10),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )
ax[0][1].errorbar(vmean, vanstd, xerr = np.std(vmean_r), yerr = np.std(vanstd_r), color = [0.0, 0.0, 0.0])
ax[0][1].annotate('Reanalysis', xy=(vmean - 0.5, vanstd + 0.1), xytext=(10, 4.5),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )

ax[1][0].errorbar(vmean, pers,   xerr = np.std(vmean_r), yerr = np.std(pers_r),   color = [0.0, 0.0, 0.0])
ax[1][0].annotate('Reanalysis', xy=(vmean - 0.5, pers + 3), xytext=(20, 35),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )

ax[1][1].errorbar(vmean, tre,    xerr = np.std(vmean_r), yerr = np.std(tre_r),    color = [0.0, 0.0, 0.0])
ax[1][1].annotate('Reanalysis', xy=(vmean - 0.5, tre + 0.05), xytext=(10, 0.2),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )


ax[0][0].grid()
ax[0][0].set_ylabel("Amplitude seasonal cycle\nvolume [10$^3$ km$^3$]", fontsize = 16)
ax[0][0].text(2.0, 18, "a", fontweight = "bold", fontsize = 16)

ax[0][1].grid()
ax[0][1].set_ylabel("Annual mean volume\nvariability [10$^3$ km$^3$]", fontsize = 16)
ax[0][1].text(2.0, 4.5, "b", fontweight = "bold", fontsize = 16)

ax[1][0].grid()
ax[1][0].set_ylabel("Persistence\nvolume [months]", fontsize = 16)
ax[1][0].text(2.0, 36, "c", fontweight = "bold", fontsize = 16)

ax[1][1].grid()
ax[1][1].set_ylabel("Trend annual mean\nvolume [10$^3$ km$^3$/y]", fontsize = 16)
ax[1][1].text(2.0, 0.35, "d", fontweight = "bold", fontsize = 16)

f = fig.colorbar(CS3, ax = ax.ravel().tolist())
years_cbar = np.arange(1930, 2090 + 1, 10)
f.set_ticks(years_cbar)
f.set_ticklabels([str(y - 10) + "-" + str(y + 10) for y in years_cbar])
fig.text(0.44, 0.02, 'Annual mean sea-ice volume [10$^3$ km$^3$]', fontsize = 18, ha='center')
plt.suptitle("CESM large ensemble (35 members) historical + RCP8.5 forcings", fontsize = 20)
plt.savefig("./fig3.png", dpi = 500)
plt.savefig("./fig3.pdf")
