#!/usr/bin/python
#
# Script to make Fig. 1 of the Process paper
# intended to introduce the context and the scientific 
# question
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


# Start the script
# ----------------
# Create the figure
fig = plt.figure(figsize = (14, 7))
gs=GridSpec(2, 4) # 2 rows, 4 columns

ax = []
ax.append(fig.add_subplot(gs[0:2, 0:2]))
ax.append(fig.add_subplot(gs[0, 2]))
ax.append(fig.add_subplot(gs[0, 3]))
ax.append(fig.add_subplot(gs[1, 2]))
ax.append(fig.add_subplot(gs[1, 3]))

# PART 1. The cycles of sea ice extent from obs
# ---------------------------------------------
# The input data. Can be retrieved directly from 
# ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/
filein = "./csv/N_seaice_extent_daily_v2.1.csv"
# Read in the data
j = 0

mydata  = list()
with open(filein, 'rb') as csvfile:
  obj = csv.reader(csvfile, delimiter = ',')
  nd = obj.line_num - 2 # nb data
  for row in obj:
    if j <= 1:
      print("ignore, header")
    else:
      mydata.append(( datetime.date(int(row[0]), int(row[1]), int(row[2])), \
                      float(row[3])))
    j = j + 1

yearmin = mydata[0][0].year
yearmax = mydata[-1][0].year
nyear   = yearmax - yearmin + 1

# Compute range on detrended anomalies
# Take the year 2015 as time reference, it has 365 data points
days_in_years = [m[0] for m in mydata if m[0].year == 2015]
ran = np.empty((365, 1))
for i in range(len(days_in_years)):
  M = days_in_years[i].month
  D = days_in_years[i].day
  a = [m[1] for m in mydata if m[0].month == M and m[0].day == D]
  # Detrend quadratically
  a_d = detrend(np.array(a), order = 1)
  # Range
  ran[i] = np.max(a_d) - np.min(a_d)

ax[0].bar([d.timetuple().tm_yday for d in days_in_years], ran, facecolor = [0.4, 0.4, 0.4], edgecolor = [0.4, 0.4, 0.4])
ax[0].text(180, 0.2, "Range (trend removed)", color = [0.9, 0.9, 0.9])
ax[0].text(5, 16.2, "a", fontweight = "bold", fontsize = 16)

# Print range of range
print("Max range: " + str(np.max(ran)))
print("Min range: " + str(np.min(ran)))
# Hack to have a color bar for each curve without making a contourf 
Z = [[0,0],[0,0]]
levels = range(yearmin, yearmax + 1, 1)
CS3 = ax[0].contourf(Z, levels, cmap = plt.cm.Reds)
c = fig.colorbar(CS3, ax = ax[0], orientation = "horizontal", pad = 0.05)
c.set_ticks(np.arange(1980, 2015 + 5, 5))
for year in np.arange(yearmin, yearmax + 1):
  data = [[m[0].timetuple().tm_yday, m[1]] for m in mydata if m[0].year == year and (m[0].month != 2 or m[0].day != 29)]
  t = [d[0] for d in data]
  col = plt.cm.Reds(int(1.0 * (year - yearmin) / nyear * 255))
  ax[0].plot(t, [d[1] for d in data], lw = 2, color = col )
ax[0].set_ylim(0, 17.0)
ax[0].set_ylabel("10$^6$ km$^2$")
ndpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # number of days per month
ax[0].xaxis.set_ticks(np.cumsum(ndpm) - ndpm[0])
ax[0].xaxis.set_ticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
for tick in ax[0].xaxis.get_majorticklabels():
    tick.set_horizontalalignment("left")
ax[0].set_title("Observed Arctic sea-ice extent, 1979-2017")
ax[0].grid()



# PART 2. The maps
# ----------------
yearb_past = 1850
yeare_past = 1880
yearb_future = 2020
yeare_future = 2050

# Load the data of SIC
ncload("./netcdfs/b.e11.B20TRC5CNBDRD.f09_g16.001.cice.h.aice_nh.185001-200512.nc", ["aice", "TLON", "TLAT"])
aice1 = aice
lon = TLON
lat = TLAT
ncload("./netcdfs/b.e11.BRCP85C5CNBDRD.f09_g16.001.cice.h.aice_nh.200601-208012.nc", ["aice"])
aice2 = aice
siconc = np.concatenate((aice1, aice2), axis = 0)

# Load the data of SIT
ncload("./netcdfs/b.e11.B20TRC5CNBDRD.f09_g16.001.cice.h.hi_nh.185001-200512.nc", ["hi"])
hi1 = hi
del hi
ncload("./netcdfs/b.e11.BRCP85C5CNBDRD.f09_g16.001.cice.h.hi_nh.200601-208012.nc", ["hi"])
hi2 = hi
sivol = np.concatenate((hi1, hi2), axis = 0)
del hi

# Create two masks where ice is always present
t1_past = (yearb_past - 1850) * 12
t2_past = (yeare_past - 1850) * 12 + 12

t1_future = (yearb_future - 1850) * 12
t2_future = (yeare_future - 1850) * 12 + 12

mask_past   = (np.min(sivol[t1_past:t2_past,     :, :], axis = 0) >= 0.0)
mask_future = (np.min(sivol[t1_future:t2_future, :, :], axis = 0) >= 0.0)

# Extract time frames of interest
owf_past = np.mean([np.max(siconc[t:t + 12, :, :], axis = 0) - \
               np.min(siconc[t:t + 12, :, :], axis = 0)   \
               for t in np.arange(t1_past, t2_past, 12)], axis = 0)

owf_past[~ mask_past] = np.nan

owf_future  = np.mean([np.max(siconc[t:t + 12, :, :], axis = 0) - \
               np.min(siconc[t:t + 12, :, :], axis = 0)   \
               for t in np.arange(t1_future, t2_future, 12)], axis = 0)
owf_future[~ mask_future] = np.nan

htp_past = np.mean([np.max(sivol[t + 12:t + 24, :, :], axis = 0) - \
               np.min(sivol[t:t + 12, :, :], axis = 0)   \
               for t in np.arange(t1_past, t2_past, 12)], axis = 0)
htp_past[~ mask_past] = np.nan

htp_future  = np.mean([np.max(sivol[t + 12:t + 24, :, :], axis = 0) - \
               np.min(sivol[t:t + 12, :, :], axis = 0)   \
               for t in np.arange(t1_future, t2_future, 12)], axis = 0)

htp_future[~ mask_future] = np.nan

summersic_past = np.mean(siconc[t1_past + 8:t2_past:12, :, :], axis = 0)
summersic_future = np.mean(siconc[t1_future + 8:t2_future:12, :, :], axis = 0)
summersic_past[summersic_past < 0.0] = np.nan
summersic_past[summersic_past > 100.0] = np.nan
summersic_future[summersic_future < 0.0] = np.nan
summersic_future[summersic_future > 100.0] = np.nan

# Plot it
# Upper left
ax[1].set_title("Open water formed\n" + str(yearb_past) + "-" + str(yeare_past) + " (model)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[1])
x, y = map(lon, lat)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, owf_past, np.arange(0, 110, 10), cmap = plt.cm.PuBu, \
                    latlon = False, extend = "neither")
map.contour(x, y, summersic_past, [15.0], colors = '#ffcccc', linewidths = 2, linestyles = 'solid')
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(0, 105, 20))
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[1].text(map(230, 65)[0], map(230, 65)[1], "b", fontweight = "bold", fontsize = 16)

# Upper right
ax[2].set_title("Open water formed\n" + str(yearb_future) + "-" + str(yeare_future) + " (model)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[2])
x, y = map(lon, lat)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, owf_future, np.arange(0, 110, 10), cmap = plt.cm.PuBu, \
                    latlon = False, extend = "neither")
map.contour(x, y, summersic_future, [15.0], colors = '#ffcccc', linewidths = 2, linestyles = 'solid')
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(0, 105, 20))
cbar.set_label("%")
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[2].text(map(230, 65)[0], map(230, 65)[1], "c", fontweight = "bold", fontsize = 16)

# Lower left
ax[3].set_title("Ice thickness seasonal\nchange " + str(yearb_past) + "-" + str(yeare_past) + " (model)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[3])
x, y = map(lon, lat)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, htp_past, np.arange(0, 2.5 + 0.25, 0.25), cmap = mycmap, \
                    latlon = False, extend = "max")
map.contour(x, y, summersic_past, [15.0], colors = '#ffcccc', linewidths = 2, linestyles = 'solid')
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(0, 2.5 + 0.5, 0.5))
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[3].text(map(230, 65)[0], map(230, 65)[1], "d", fontweight = "bold", fontsize = 16)

# Lower right
ax[4].set_title("Ice thickness seasonal\nchange " + str(yearb_future) + "-" + str(yeare_future) + " (model)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[4])
x, y = map(lon, lat)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, htp_future, np.arange(0, 2.5 + 0.25, 0.25), cmap = mycmap, \
                    latlon = False, extend = "max")
map.contour(x, y, summersic_future, [15.0], colors = '#ffcccc', linewidths = 2, linestyles = 'solid')
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(0, 2.5 + 0.5, 0.5))
cbar.set_label("m")
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[4].text(map(230, 65)[0], map(230, 65)[1], "e", fontweight = "bold", fontsize = 16)


# SAVE FIGURE
# -----------
plt.savefig("./fig1.png", dpi = 500)
plt.savefig("./fig1.pdf")



