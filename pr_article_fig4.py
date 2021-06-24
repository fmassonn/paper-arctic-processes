#!/usr/bin/python
from scipy.signal import gaussian
from scipy.ndimage import filters
# Load namelist
exec(open("./config.py").read())

# Start the script
domain = "Arctic_n80" # Domain under consideration (north of X)
dirin = "./netcdfs/"   # where the pre-processed NetCDF data of volume and area are

# Years spanning the period investigated (time vector)
yearb = 1850
yeare = 2300

# Years defining period to calculate trends
yearb_f = 2020
yeare_f = 2050

# Years definind the refence period for mean state
yearb_r = 2000
yeare_r = 2010

years = np.arange(yearb, yeare + 1)
n_year = yeare - yearb + 1
time = np.arange(yearb, yeare + 1, 1.0 / 12)


plt.figure("fig1", figsize = (5, 3))
plt.figure("fig2", figsize = (4, 4))
plt.figure("fig3", figsize = (4, 4))

# PART 1 - Load observational data
# --------------------------------
plt.figure("fig1")
# Will have the 2 obs
obs = np.empty((len(time), 3))
obs[:] = np.nan

# Icesat
col = [0.8, 0.2, 0.2]
ma = '*'
ncload(dirin + "icesat_200301-200812.nc", "sivol_" + domain)
t1 = (2003 - yearb) * 12
t2 = (2008 - yearb) * 12 + 12
obs[t1:t2, 0] = eval("sivol_" + domain)
plt.scatter(time, obs[:, 0], marker = ma, color = col, zorder = 20)
plt.scatter(2001, 19.0, marker = ma, color = col)
plt.text(2001, 18.5, "  Estimate IceSat", color = col)

# CryoSat2
col = [1.0, 0.5, 0.5]
ma = 'v'
ncload(dirin + "cryosat2_201001-201712.nc", "sivol_" + domain)
t1 = (2010 - yearb) * 12
t2 = (2017 - yearb) * 12 + 12
obs[t1:t2, 1] = eval("sivol_" + domain)
plt.scatter(time, obs[:, 1], marker = ma, color = col, zorder = 15)
plt.scatter(2008, 19.0, marker = ma, color = col)
plt.text(2008, 18.5, "  Estimate CryoSat2", color = col)

# ITRP
col = [1.0, 0.7, 0.7]
ma = 'o'
ncload(dirin + "itrp_200001-201212.nc", "sivol_" + domain)
t1 = (2000 - yearb) * 12
t2 = (2012 - yearb) * 12 + 12
obs[t1:t2, 2] = eval("sivol_" + domain)
plt.scatter(time, obs[:, 2], marker = ma,color =  col, zorder = 10)
plt.scatter(2008, 17.5, marker = ma, color = col)
plt.text(2008, 17.0, "  Estimate ITRP", color = col)



# Lower and upper bounds
obs_lb, obs_ub = np.nanmin(obs, axis = 1), np.nanmax(obs, axis = 1)

# PART 2 - Load and analyse CMIP5 data
# ------------------------------------

# Read CMIP5 namelist (in another file since this namelist is used by different scripts)
# Read historical namelist
exec(open("./cmip5_namelist.py").read())
info_historical = info
execfile("./cmip5_namelist_RCP85.py")
info_rcp85 = info

n_models_historical = len(info_historical)
n_models_rcp85      = len(info_rcp85)

data = list()
reject = list()
dV = list()

np.random.seed(1) # to fix the colors of CMIP5
colors = [np.random.random(3) for j in range(n_models_historical)]


n_models_match = 0

for j_models_historical in range(n_models_historical):
  # 1. Get information about the model
  model_historical = info_historical[j_models_historical][0]
  yb_historical    = info_historical[j_models_historical][5]
  ye_historical    = info_historical[j_models_historical][6]
  n_members_historical = len(info_historical[j_models_historical][7])
  # Let us find if there is a match with the RCP8.5 namelist
  for j_models_rcp85 in range(n_models_rcp85):
    model_rcp85 = info_rcp85[j_models_rcp85][0]
    n_members_rcp85 = len(info_rcp85[j_models_rcp85][7])
    yb_rcp85 = info_rcp85[j_models_rcp85][5]
    ye_rcp85 = info_rcp85[j_models_rcp85][6]


    if model_rcp85 == model_historical:

      # Loop through members of the historical simu and attempt to find match
      for j_members_historical in np.arange(n_members_historical):
        member_historical = info_historical[j_models_historical][7][j_members_historical]

        for j_members_rcp85 in np.arange(n_members_rcp85):
          member_rcp85 = info_rcp85[j_models_rcp85][7][j_members_rcp85]

          if member_rcp85 == member_historical:
            print(str(model_historical) + " r" + str(member_historical) + "i1p1: MATCH!!")
            n_models_match += 1
            model_match = True

            # Load the historical simulation
            n_models_match += 1
            model_match = True

            # Load the historical simulation
      # Loop through members of the historical simu and attempt to find match
      for j_members_historical in np.arange(n_members_historical):
        member_historical = info_historical[j_models_historical][7][j_members_historical]

        for j_members_rcp85 in np.arange(n_members_rcp85):
          member_rcp85 = info_rcp85[j_models_rcp85][7][j_members_rcp85]

          if member_rcp85 == member_historical:
            print(str(model_historical) + " r" + str(member_historical) + "i1p1: MATCH!!")
 
            # Load the historical simulation
            file_historical = dirin + "seaice_" + model_historical + "_historical_r" + str(member_historical) + "i1p1_" + str(yb_historical) + "01-" + str(ye_historical) + "12.nc"

            ncload(file_historical, "sivol_" + domain)

            tmp = eval("sivol_" + domain)
            
            series_record = np.empty((n_year * 12)) # Vector with same length as time axis
            series_record[:] = np.nan
 
            series_record[(yb_historical - yearb) * 12 : (ye_historical - yearb) * 12 + 12] = tmp

            del tmp

            # Load the RCP8.5 simulation and override if necessary
            file_rcp85 = dirin + "seaice_" + model_rcp85 + "_rcp85_r" + str(member_rcp85) + "i1p1_" + str(yb_rcp85) + "01-" + str(ye_rcp85) + "12.nc"
            ncload(file_rcp85, "sivol_" + domain)

            tmp = eval("sivol_" + domain)

            series_record[(yb_rcp85 - yearb) * 12 : (ye_rcp85 - yearb) * 12 + 12] = tmp

            # Append the result
            data.append(series_record)

            # Make analysis
            plt.figure("fig1")
            plt.plot(time, series_record, color = colors[j_models_historical])  


            # Make test of compatibility
            if (np.nanmin(series_record - obs_ub) > 0.0 and np.nanmax(series_record - obs_ub) > 0.0) \
               or \
               (np.nanmin(series_record - obs_lb) < 0.0 and np.nanmax(series_record - obs_lb) < 0.0):

              reject.append(True)
              color = [0.6, 0.6, 0.6]
              lw = 0.5
              zorder = 0
            else:
              color = colors[j_models_historical]
              lw = 1.0
              zorder = 1000.0
              reject.append(False)
              print("")
              print(" >>>>> " + model_historical)
              print("")
            
            # Compute annual means
            anmean = [np.mean(series_record[t:t + 12]) for t in np.arange(0, len(series_record), 12)]
           
            plt.figure("fig3")
            #b = gaussian(5, 10)

            #ga = filters.convolve1d(anmean, b/b.sum())
            ga = anmean
            plt.plot(years, anmean ,color = color, lw = lw, zorder = zorder)


            # Compute future loss
            loss = np.polyfit(np.arange(yearb_f, yeare_f + 1), anmean[(yearb_f - yearb) : (yeare_f - yearb + 1)], 1)[0] * (yeare_f - yearb_f + 1)
            dV.append(loss)

            # Scatter plot trends
            mean_state = np.nanmean(series_record[(yearb_r - yearb) * 12 : (yeare_r - yearb) * 12 + 12])
            plt.figure("fig2")
            plt.scatter(mean_state, loss, 10, color = colors[j_models_historical])

            del series_record

dV_sub = [dV[i] for i in range(len(reject)) if reject[i] == False]

plt.figure("fig1")
plt.grid()
plt.xlim(2000, 2020)
plt.ylim(0.0, 20.0)
plt.ylabel("Sea-ice volume >80N [10$^3$ km$^3$]")
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()
plt.savefig("./figED16.png", dpi = 300)

plt.figure("fig2")
plt.xlim(0, 18)
plt.xlabel(str(yearb_r) + "-" + str(yeare_r) + " annual mean\nsea-ice volume >80N [10$^3$ km$^3$]")
plt.ylabel(str(yearb_f) + "-" + str(yeare_f) + " annual mean\nsea-ice volume change >80N [10$^3$ km$^3$]")
plt.grid()
plt.title(str(len(dV)) + " CMIP5 simulations\n(historical + RCP8.5 scenarios)")
plt.tight_layout()
plt.savefig("./figED15.png", dpi = 300)
plt.clf()


plt.figure("fig3")
plt.grid()
plt.ylim(0.0, 20.0)
plt.xlim(2000, 2100)
plt.title(str(len(dV)) + " CMIP5 simulations\n(historical + RCP8.5 scenarios)")
plt.text(2010, 18.5, str(yearb_f) + "-" + str(yeare_f) + " ice-volume loss (>80N):")
plt.text(2010, 17, "ALL: " + str(np.round(np.mean(dV), 2)) + " +/- " + str(np.round(np.std(dV), 2)) + " x 10$^3$ km$^3$", color = [0.6, 0.6, 0.6])
plt.text(2010, 16, "SUB: " + str(np.round(np.mean(dV_sub), 2)) + " +/- " + str(np.round(np.std(dV_sub), 2)) + " x 10$^3$ km$^3$", color = [0.2, 0.7, 0.2])
plt.ylabel("Annual mean sea-ice volume [10$^3$ km$^3$]")
plt.tight_layout()
plt.savefig("./fig4.png", dpi = 300)
plt.savefig("./fig4.pdf")
plt.clf()
