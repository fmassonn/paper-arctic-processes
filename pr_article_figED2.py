#!/usr/bin/python
#
# Script to make Fig. S2 of the process paper
# intended to show the winter/summer year-to-year
# changes in sea ice extent per season
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


# Start the script
# ----------------
# Create the figure
fig = plt.figure(figsize = (5, 3))

# Load data
ncload("netcdfs/siextentn_r1i1p1_mon_197801-201712.nc", "siextentn")

# Get March and September time series
march_sie = siextentn[2::12]
sept_sie  = siextentn[8::12]

years = np.arange(1978, 2017 + 1)


plt.plot(years, march_sie, lw = 2, color = [0.5, 0.5, 1.0])
plt.plot(years, sept_sie, lw = 2, color = [1, 0.5, 0.5])

plt.text(1980, 17, "March extent", color = [0.5, 0.5, 1.0])
plt.text(1980, 10, "September extent", color = [1.0, 0.5, 0.5])

w2s = sept_sie - march_sie
s2w = march_sie[1:] - sept_sie[:-1]

for year in years[:-1]:
  # Winter to summer loss
  shade = (np.abs(w2s[year - 1978]) - np.min(np.abs(w2s))) / (np.max(np.abs(w2s)) - np.min(np.abs(w2s)))
  plt.bar(year - 0.4, w2s[year - 1978], width = 0.8, color = [1.0 - shade * 0.3, 1.0 - shade, 1.0 - shade], edgecolor = [0.5, 0.5, 0.5], lw = 0.5)
  # Summer to winter gain
  shade = (np.abs(s2w[year - 1978]) - np.min(np.abs(s2w))) / (np.max(np.abs(s2w)) - np.min(np.abs(s2w))) 
  plt.bar(year + 0.5 - 0.4, s2w[year - 1978], width = 0.8, color = [1.0 - shade, 1.0 - shade, 1.0 - shade * 0.3], edgecolor = [0.5, 0.5, 0.5], lw = 0.5)

plt.text(2018, 4, "September\nto March gain", color = [0.5, 0.5, 1.0], fontsize = 10)
plt.text(2018, - 7, "March to\nSeptember loss", color = [1.0, 0.5, 0.5], fontsize = 10)
plt.ylabel("10$^6$ km$^2$")
plt.title("Seasonal budget of sea-ice extent changes")
plt.xlim(1975, 2032)
# SAVE FIGURE
# -----------
plt.grid()
plt.tight_layout()
plt.savefig("./figED2.png", dpi = 500)



