#!/usr/bin/python
#
# Script to make ED Fig. 1 of the process paper
# intended to show the winter/summer ratio of 
# temperature trends
#
# Author: F. Massonnet, September 2017
#         francois.massonnet@uclouvain.be

# Load standard modules, fonts, etc.
execfile("./config.py") 


# Start the script
# ----------------
# Create the figure
fig = plt.figure(figsize = (13.5, 4.5))
gs=GridSpec(1, 3) # 2 rows, 4 columns

ax = []
ax.append(fig.add_subplot(gs[0]))
ax.append(fig.add_subplot(gs[1]))
ax.append(fig.add_subplot(gs[2]))


# Load data
ncload("./netcdfs/tas_r1i1p1_mon_197901-201612_reg0.75.nc", ["tas", "longitude", "latitude"])
longitude, latitude = np.meshgrid(longitude, latitude) # Expand to 2-D grid
nt, ny, nx = tas.shape

# Compute JFM averages
tas_jfm = np.array([np.mean(tas[t + 0:t + 3, :, :], axis = 0) for t in np.arange(0, nt, 12)])
# Compute JAS averages
tas_jas = np.array([np.mean(tas[t + 6:t + 9, :, :], axis = 0) for t in np.arange(0, nt, 12)])

nyear = tas_jfm.shape[0]
# Compute JFM trends
tas_jfm_change = np.reshape(np.polyfit(range(nyear), np.reshape(tas_jfm, (nyear, ny * nx), order = 'C'), 1)[0], (ny, nx), order = 'C') * nyear
tas_jas_change = np.reshape(np.polyfit(range(nyear), np.reshape(tas_jas, (nyear, ny * nx), order = 'C'), 1)[0], (ny, nx), order = 'C') * nyear

# Make maps
# Add one column to longitude to avoid white stripe
longitude = np.concatenate(     (longitude, longitude[:, :1]          ), axis = 1)
latitude = np.concatenate(      (latitude, latitude[:, :1]            ), axis = 1)
tas_jfm_change = np.concatenate((tas_jfm_change, tas_jfm_change[:, :1]), axis = 1)
tas_jas_change = np.concatenate((tas_jas_change, tas_jas_change[:, :1]), axis = 1)

ax[0].set_title("Winter (JFM) near-surface temperature\nchange (1979-2016, reanalysis)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[0])
x, y = map(longitude, latitude)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, tas_jfm_change, np.arange(-10, 10, 1), cmap = plt.cm.RdBu_r, \
                    latlon = False, extend = "both")
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(-10, 10 + 2, 2))
cbar.set_label("$^\circ$C")
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[0].text(map(230, 65)[0], map(230, 65)[1], "a", fontweight = "bold", fontsize = 16)


ax[1].set_title("Summer (JAS) near-surface temperature\nchange (1979-2016, reanalysis)")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[1])
x, y = map(longitude, latitude)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, tas_jas_change, np.arange(-10, 10, 1), cmap = plt.cm.RdBu_r, \
                    latlon = False, extend = "both")
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(-10, 10 + 2, 2))
cbar.set_label("$^\circ$C")
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[1].text(map(230, 65)[0], map(230, 65)[1], "b", fontweight = "bold", fontsize = 16)

ax[2].set_title("Ratio")
map = Basemap(projection = "npstere", boundinglat = 70.0, \
                lon_0 = 0.0, resolution = 'l', ax = ax[2])
x, y = map(longitude, latitude)
map.drawcoastlines(linewidth = 0.25)
map.fillcontinents(color = 'grey', lake_color = 'w')
cs = map.contourf(x, y, np.abs(tas_jfm_change / tas_jas_change), np.arange(0, 5 + 0.5, 0.5), cmap = plt.cm.YlOrBr, \
                    latlon = False, extend = "both")
cbar = map.colorbar(cs, location = 'right', pad = "5%")
cbar.set_ticks(np.arange(0, 6, 1))
cbar.cmap.set_over((plt.cm.YlOrBr(255)[:3]))
cbar.set_label("[-]")
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,10))
ax[2].text(map(230, 65)[0], map(230, 65)[1], "c", fontweight = "bold", fontsize = 16)


# SAVE FIGURE
# -----------
plt.tight_layout()
plt.savefig("./figED1.png", dpi = 500)



