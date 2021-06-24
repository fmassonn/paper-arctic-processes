#!/usr/bin/python
#
# Script to make ED Fig. 10 of the Process paper
# intended to explain the general methodology
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


# Start the script
domain = "Arctic_n80" # Domain under consideration (north of X)
dirin = "./netcdfs/"   # where the pre-processed NetCDF data of volume and area are

# Years defining the period to showcase (must be within the availability of data)
yearb = 1955
yeare = 2004

# Read CMIP5 namelist (in another file since this namelist is used by different scripts)
execfile("./cmip5_namelist.py")

# ID of the models to display
id_model = 0


# Start the figure
# Create a subplot with two rows (one for volume, one for area)
# and three columns with uneven size: first column to plot time series
# and second and third to show scatter plots from two models

fig=plt.figure(figsize = (9, 7))
gs=GridSpec(2, 3) # 2 rows, 3 columns

ax1 = fig.add_subplot(gs[0, 0:2]) # First row, first column
ax2 = fig.add_subplot(gs[0, 2]  ) # First row, second column
ax4 = fig.add_subplot(gs[1, 0:2]) # Second row, first column
ax5 = fig.add_subplot(gs[1, 2]  ) # Second row, second column

# Loop over models and plot data from the models's first member
model = info[id_model][0]
member = "r1i1p1"
yb = info[id_model][5]
ye = info[id_model][6]

ncload(dirin + "/seaice_" + model + "_historical_" + member + "_" + str(yb) + "01-" + str(ye) + "12.nc", ["sivol_" + domain, "siarea_" + domain])

volume = eval("sivol_" + domain)
area   = eval("siarea_" + domain)

# Time axis
base = datetime.datetime(yb, 1, 1)
time = [base + datetime.timedelta(days = jt * 365.0 / 12.0) for jt in range(len(volume))]

# Subset to period of interest
t1 = (yearb - yb) * 12 
t2 = (yeare - yb) * 12 + 12

 # Time of minimum and maximum of sea ice volume
tmin = [jt + np.argmin(volume[jt:jt + 12]) for jt in np.arange(t1, t2 - 12, 12)]
tmax = [jt + np.argmax(volume[jt:jt + 12]) for jt in np.arange(t1, t2, 12)]

# Time series of interest
vmin = np.array([volume[t] for t in tmin])
vmax = np.array([volume[t] for t in tmax])
dv_growth   = vmax[1:] - vmin

amin = np.array([area[t] for t in tmin])
amax = np.array([area[t] for t in tmax])
da_melt   = amin - amax[:-1] 
dv_melt   = vmin - vmax[:-1]

# Detrended (order 1) values
dv_growth_d = detrend(dv_growth, order = 1)
vmin_d      = detrend(vmin     , order = 1)
da_melt_d   = detrend(da_melt  , order = 1)
dv_melt_d   = detrend(dv_melt  , order = 1)

# Plot time series of volume
# --------------------------
# Colors found here: http://www.colorhexa.com/
ax1.plot(time[t1:t2], volume[t1:t2], lw = 2, color = "#537d8d")
# Locate minima and maxima
# Plot minima and maxima
[ax1.scatter(time[t], volume[t], 50, color = "#203036") for t in tmin]
[ax1.scatter(time[t], volume[t], 50, color = "#aec6cf") for t in tmax]
ax1.set_ylabel("10$^3$ km$^3$")
ax1.set_xlim(datetime.datetime(1970, 1, 1, 0, 0), datetime.datetime(1980, 1, 1, 0, 0))
ax1.set_ylim(0, 11)
ax1.set_title("Sea-ice volume north of 80N (" + model + ")")
ax1.text(datetime.datetime(1970, 2, 1, 0, 0), 10.2, "a", fontweight = "bold", fontsize = 16)
ax1.grid()


# Repeat with area
# --------------------------
ax4.plot(time[t1:t2], area[t1:t2], lw = 2, color = "#ff8243")
# Plot minima and maxima
[ax4.scatter(time[t], area[t], 50, color = "#f65200") for t in tmin]
[ax4.scatter(time[t], area[t], 50, color = "#ffb590") for t in tmax]
ax4.set_ylabel("10$^6$ km$^2$")
ax4.set_xlim(datetime.datetime(1970, 1, 1, 0, 0), datetime.datetime(1980, 1, 1, 0, 0))
ax4.set_title("Sea-ice area north of 80N (" + model + ")")
ax4.set_ylim(0, 4)
ax4.text(datetime.datetime(1970, 2, 1, 0, 0), 3.7, "c", fontweight = "bold", fontsize = 16)
ax4.grid()

# Scatter plot the detrended values and mark regression line
ax2.scatter(vmin_d, dv_growth_d, 10, color = "#537d8d", marker = ',')
ax2.plot((-1e9, 1e9), (0, 0), lw = 1, color = 'k')
ax2.plot((0, 0), (-1e9, 1e9), lw = 1, color = 'k')
fit = np.polyfit(vmin_d, dv_growth_d, 1)
xfit = np.arange(-10, 10)
ax2.plot(xfit, fit[0] * xfit + fit[1], color = "#206060", linestyle = '-', linewidth = 2)
ax2.set_xlabel("Volume at minimum\n(anomalies) [10$^3$ km$^3$]")
ax2.set_ylabel("Wintertime volume range\n(anomalies) [10$^3$ km$^3$]")
ax2.set_xlim(-2.5, 2.5)
ax2.set_ylim(-2.5, 2.5)
ax2.set_title("Evaluation of the\nIFE (" + str(yearb) + "-" + str(yeare) + ")")
ax2.text(-2.0, -2.0, "y = " + str(np.round(fit[0], 2)) + " x", color = "#206060")
ax2.text(-2.4, 2.1, "b", fontweight = "bold", fontsize = 16)

ax2.grid()

ax5.scatter(dv_melt_d, da_melt_d, 10, color = "#ff8243", marker = ',')
ax5.plot((-1e9, 1e9), (0, 0), lw = 1, color = 'k')
ax5.plot((0, 0), (-1e9, 1e9), lw = 1, color = 'k')
fit = np.polyfit(dv_melt_d, da_melt_d, 1)
xfit = np.arange(-10, 10)
ax5.plot(xfit, fit[0] * xfit + fit[1], color = "#ff6f00", linestyle = '-', linewidth = 2)
ax5.set_xlabel("Summertime volume range\n(anomalies) [10$^3$ km$^3$]")
ax5.set_ylabel("Summertime area range\n(anomalies) [10$^6$ km$^2$]")
ax5.set_xlim(-2.5, 2.5)
ax5.set_ylim(-2.5, 2.5)
ax5.set_title("Evaluation of the\nOWFE (" + str(yearb) + "-" + str(yeare) + ")")
ax5.text(-2.0, -2.0, "y = " + str(np.round(fit[0], 2)) + " x", color = "#ff6f00")
ax5.text(-2.4, 2.1, "d", fontweight = "bold", fontsize = 16)
ax5.grid()


plt.tight_layout()
fig.savefig('figED10.png', dpi = 500)
