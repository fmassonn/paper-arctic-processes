#!/usr/bin/python
#
# Script to make ED Fig. 3 of the process paper
# intended to show the seasonality of CESM
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 

dirin = "./netcdfs/"
# Start the script
# ----------------
# Create the figure
fig = plt.figure(figsize = (6, 6))
gs=GridSpec(1, 1) # 2 rows, 4 columns

ax = []
ax.append(fig.add_subplot(gs[0, 0]))

# PART 1. The cycles of sea ice extent from model
# -----------------------------------------------
ncload(dirin + "b.e11.B20TRC5CNBDRD.f09_g16.001.185001-200512.nc", ["siextent_Arctic"])
extent1 = siextent_Arctic
ncload(dirin + "b.e11.BRCP85C5CNBDRD.f09_g16.001.200601-210012.nc", ["siextent_Arctic"])
extent2 = siextent_Arctic

yearmin = 1979
yearmax = 2017

extent = np.concatenate((extent1, extent2), axis = 0)
extent = extent[(yearmin - 1850) * 12:(yearmax - 1850) * 12 + 12]
nyear = yearmax - yearmin + 1

# Compute range on detrended anomalies
# Take the year 2015 as time reference, it has 365 data points
ran = np.empty((12, 1))

for m in range(12): # for each month
  a = extent[m::12]
  a_d = detrend(a, order = 1)
  ran[m] = np.max(a_d) - np.min(a_d)

ax[0].bar(np.arange(12), ran, width = 1.0, facecolor = [0.4, 0.4, 0.4], edgecolor = [0.4, 0.4, 0.4])
ax[0].text(6, 0.2, "Range (trend removed)", color = [0.9, 0.9, 0.9])

# Print range of range
#print("Max range: " + str(np.max(ran)))
#print("Min range: " + str(np.min(ran)))
# Hack to have a color bar for each curve without making a contourf 
Z = [[0,0],[0,0]]
levels = range(yearmin, yearmax + 1, 1)
CS3 = ax[0].contourf(Z, levels, cmap = plt.cm.Reds)
c = fig.colorbar(CS3, ax = ax[0], orientation = "horizontal", pad = 0.05)
c.set_ticks(np.arange(1980, 2015 + 5, 5))
for year in np.arange(yearmin, yearmax + 1):
  data = extent[(year - yearmin) * 12:(year - yearmin) * 12 + 12]
  t = np.arange(0.5, 12.5)
  col = plt.cm.Reds(int(1.0 * (year - yearmin) / nyear * 255))
  ax[0].plot(t, data, lw = 2, color = col )
ax[0].set_ylim(0, 17.0)
ax[0].set_ylabel("10$^6$ km$^2$")
ax[0].xaxis.set_ticks(np.arange(12))
ax[0].xaxis.set_ticklabels(["  Jan", "  Feb", "  Mar", "  Apr", "  May", "  Jun", "  Jul", "  Aug", "  Sep", "  Oct", "  Nov", "  Dec"])
for tick in ax[0].xaxis.get_majorticklabels():
    tick.set_horizontalalignment("left")
ax[0].set_title("CESM-LE Arctic sea-ice extent, 1979-2017 (member 001)")
ax[0].grid()


# SAVE FIGURE
# -----------
plt.tight_layout()
plt.savefig("./figED3.png", dpi = 500)



