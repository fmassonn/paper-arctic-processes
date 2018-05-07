#!/usr/bin/python
#
# Script to make Fig. 2 of the Process paper
# intended to highlight the strong dependence
# of the process factors to the mean state
#
# Also makes ED Table 1 (legend) and ED Figure 5 (CESM processes)

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
fig=plt.figure(figsize = (8, 6))
gs=GridSpec(2, 2) # 2 rows, 2 columns

ax = []
ax.append(fig.add_subplot(gs[0, 0]))
ax.append(fig.add_subplot(gs[1, 0]))
ax.append(fig.add_subplot(gs[0, 1]))
ax.append(fig.add_subplot(gs[1, 1]))


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

    data[j_models].append([vmean, nfb, pfb])


# Scatter plots
# -------------
comp = info[j_models][9] - 1 # Complexity index
fit = list()
std = list()

for j_models in range(n_models):
  members = info[j_models][7]
  n_members = len(members)

  # Plot mean of diagnostics
  x  = np.mean([i[0] for i in data[j_models]])
  y1 = np.mean([i[1] for i in data[j_models]])
  y2 = np.mean([i[2] for i in data[j_models]])
 
  ax[0].scatter(x, y1, 20 * n_members, color = colors[j_models], alpha = 0.5)
  ax[0].text(   x, y1, str(j_models + 1).zfill(2), color = [0.5 * c for c in colors[j_models]], fontsize = 4, ha = 'center', va = 'center')
  ax[1].scatter(x, y2, 20 * n_members, color = colors[j_models], alpha = 0.5)
  ax[1].text(   x, y2, str(j_models + 1).zfill(2), color = [0.5 * c for c in colors[j_models]], fontsize = 4, ha = 'center', va = 'center')
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

ax[1].set_xlabel("Annual mean ice volume [10$^3$ km$^3$]")
ax[1].set_ylabel("OWFE [m$^{-1}$]")
ax[1].plot((-1e9, 1e9), (0, 0), 'k-', lw = 3)
ax[1].plot((0, 0), (-1e9, 1e9), 'k-', lw = 3)
ax[1].set_xlim(-1.0, 25.0)
ax[1].set_ylim(-0.2, 1.0 )
ax[1].text(0.7, 0.9, "c", fontweight = "bold", fontsize = 16)
ax[1].grid()
  
# Compute 1/x fits
# All members are taken as one data point to weigh the data towards models with more
# information (otherwise the model means would be given the same weight, which is not correct)

x =  np.array([np.mean([data[j_models][j][0] for j in range(len(data[j_models]))]) for j_models in range(n_models)])
y1 = np.array([np.mean([data[j_models][j][1] for j in range(len(data[j_models]))]) for j_models in range(n_models)])
y2 = np.array([np.mean([data[j_models][j][2] for j in range(len(data[j_models]))]) for j_models in range(n_models)])
u = 1.0 / x # The data to be fitted with a linear model
fit.append(np.polyfit(u, y1, 1))   # fit are the fit parameters
residuals =  y1 - (fit[0][0] * u + fit[0][1])
std.append(np.std(residuals))
fit.append(np.polyfit(u, y2, 1))   # fit are the fit parameters
residuals =  y2 - (fit[1][0] * u + fit[1][1])
std.append(np.std(residuals))

# Print correlations with inverse of thickness
r1 = np.corrcoef(1.0 / x, y1)[0, 1]
r2 = np.corrcoef(1.0 / x, y2)[0, 1]
print("IFE: R^2 = " + str(np.corrcoef(1.0 / x, y1)[0, 1] ** 2))
print("OWFE: R^2 = " + str(np.corrcoef(1.0 / x, y2)[0, 1] ** 2))

N = len(x)
for r in [r1, r2]:
  print("R^2 = " + str(r ** 2))
  tstat = r / np.sqrt((1 - r ** 2) / (N - 2))  # The t-statistic.
                                               # Under the null hypothesis of no correlation,
                                               # tstat follows a student's law with  N - 2 dof.
  pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), N - 2)

  print("P-value (one-sided): " + str(pval))

# Add obs/reanalysis estimates
# Add OBS estimates
ncload(dirin + "icesat_200301-200812.nc", ["sivol_" + domain])
ncload(dirin + "osisaf_200301-200812.nc", ["siarea_" + domain])

volume = eval("sivol_" + domain)
area   = eval("siarea_" + domain)

nf, stats, _ = negative_seaice_feedback(volume, period = 12)
# Mean state and 1 std
vm = np.nanmean(volume)
vm_std = np.nanstd(volume)

ax[0].scatter(vm, nf, 40, marker = "x", color = [0.5, 0.2, 0.0])
ax[0].errorbar(vm, nf, xerr = vm_std, yerr = stats[2], color = [0.5, 0.2, 0.0], zorder = 50)
ax[0].text(7, -0.85, "OBS", color = [0.5, 0.2, 0.0])

pf, stats, _ = positive_seaice_feedback(volume, area, period = 12)
ax[1].scatter(vm, pf, 40, marker = "x", color = [0.5, 0.2, 0.0])
ax[1].errorbar(vm, pf, xerr = vm_std, yerr = stats[2], color = [0.5, 0.2, 0.0], zorder = 50)
ax[1].text(4.75, 0.02, "OBS", color = [0.5, 0.2, 0.0])

del volume, area, nf, stats, pf, vm, vm_std

# Add PIOMAS estimates
ncload(dirin + "piomas_197901-201512.nc", ["sivol_" + domain, "siarea_" + domain])
volume = eval("sivol_" + domain)
area = eval("siarea_" + domain)

vm = np.nanmean(volume)
vanmean = np.array([np.mean(volume[t: t + 12]) for t in np.arange(0, len(volume), 12)])
vm_std = np.nanstd(vanmean)

nf, stats, _ = negative_seaice_feedback(volume, period = 12)
ax[0].scatter(vm, nf, 40, marker = "*", color = [0.0, 0.0, 0.0])
ax[0].errorbar(vm, nf, xerr = vm_std, yerr = stats[2], color = [0.0, 0.0, 0.0], zorder = 50)
del stats
ax[0].annotate('Reanalysis', xy=(vm + 0.5, nf - 0.05), xytext=(10, -0.7),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )

pf, stats, _ = positive_seaice_feedback(volume, area, period = 12)
ax[1].scatter(vm, pf, 40, marker = "*", color = [0.0, 0.0, 0.0])
ax[1].errorbar(vm, pf, xerr = vm_std, yerr = stats[2], color = [0.0, 0.0, 0.0], zorder = 50)
del stats
ax[1].annotate('Reanalysis', xy=(vm + 0.5, pf + 0.05), xytext=(12, 0.4),
            arrowprops = dict(facecolor = [0.0, 0.0, 0.0], shrink=0.05, width = 0.1, headwidth = 2, headlength = 2),
            )


del volume, area


# PART 2 - Load and analyse toy model data
# ----------------------------------------
exec(open("./sim1D_namelist.py").read())
n_exps = len(exps)

data = list()

for j_exps in range(n_exps):
  data.append(list())
  yearb = exps[j_exps][1]
  yeare = exps[j_exps][2]
  ncload(dirin + "sim1D_" + exps[j_exps][0] + "_18500101-19491231.nc", ["sivol_" + domain, "siarea_" + domain])

  # Subset to last 50 years
  volume = sivol_Arctic_n80[-365 * 50:] 
  area   = siarea_Arctic_n80[-365 * 50:] 

  # Negative fb
  nf, _, _ = negative_seaice_feedback(volume, period = 365)

  # Positive fb
  pf, _, _ = positive_seaice_feedback(volume, area, period = 365)

  # Annual-mean volume
  data[j_exps].append(np.mean(volume))
  data[j_exps].append(nf             )
  data[j_exps].append(pf             )

  del volume, area, nf, pf

# Scatter plot the result
for d, e in zip(data, exps):
  ax[2].scatter(d[0], d[1], 50, marker = e[5], color = e[4])
  ax[2].plot((-1e9, 1e9), (0, 0), 'k-', lw = 3)
  ax[2].plot((0, 0), (-1e9, 1e9), 'k-', lw = 3)
  ax[2].set_xlim(-1.0, 25.0)
  ax[2].set_ylim(-1.0, 0.2)
  ax[2].grid()
  ax[2].set_title("1-D thermodynamic model")

for d, e in zip(data, exps):
  ax[3].scatter(d[0], d[2], 50, marker = e[5], color = e[4])
  ax[3].set_xlabel("Annual mean ice volume [10$^3$ km$^3$]")
  ax[3].plot((-1e9, 1e9), (0, 0), 'k-', lw = 3)
  ax[3].plot((0, 0), (-1e9, 1e9), 'k-', lw = 3)
  ax[3].set_xlim(-1.0, 25.0)
  ax[3].set_ylim(-0.2, 1.0)
  ax[3].grid()

ax[3].scatter(10, 0.9, 50, marker = exps[0][5], color = exps[0][4])
ax[3].text(   11, 0.9, "1-D model reference", va = 'center')
ax[3].scatter(10, 0.8, 50, marker = exps[1][5], color = exps[1][4])
ax[3].text(   11, 0.8, "ice conductivity varied", color = exps[1][4], va = 'center')
ax[3].scatter(10, 0.7, 50, marker = exps[11][5], color = exps[11][4])
ax[3].text(   11, 0.7, "ice albedo varied", color = exps[11][4], va = 'center')
ax[3].scatter(10, 0.6, 50, marker = exps[21][5], color = exps[21][4])
ax[3].text(   11, 0.6, "forcing varied", color = exps[21][4], va = 'center')

# Plot CMIP5 envelopes
xx  = np.arange(1.0, 30.0, 0.05)
uu  = 1.0 / xx
yy1 = (fit[0][0] * uu + fit[0][1])
ax[2].fill_between(xx, yy1 - 2 * std[0], yy1 + 2 * std[0], linewidth = 0, edgecolor = None, facecolor = [0.8, 0.8, 0.8], zorder = 0)
ax[2].fill_between(xx, yy1 -     std[0], yy1 +     std[0], linewidth = 0, edgecolor = None, facecolor = [0.7, 0.7, 0.7], zorder = 0)
ax[2].text(15, -0.65, "CMIP5", color = [0.7, 0.7, 0.7])
ax[2].text(0.7, 0.1, "b", fontweight = "bold", fontsize = 16)


yy2 = (fit[1][0] * uu + fit[1][1])
ax[3].fill_between(xx, yy2 - 2 * std[1], yy2 + 2 * std[1], linewidth = 0, edgecolor = None, facecolor = [0.8, 0.8, 0.8], zorder = 0)
ax[3].fill_between(xx, yy2 -     std[1], yy2 +     std[1], linewidth = 0, edgecolor = None, facecolor = [0.7, 0.7, 0.7], zorder = 0)
ax[3].text(15, 0.4, "CMIP5", color = [0.7, 0.7, 0.7])
ax[3].text(0.7, 0.9, "d", fontweight = "bold", fontsize = 16)

plt.tight_layout()
plt.savefig("./fig2.png", dpi = 500)
plt.savefig("./fig2.pdf")

# Print model names in colors
fig = plt.figure(figsize = (4.7, 8))
cplx_lab = ["Very simple", "Simple", "Intermediate", "Complex"]
plt.text(0, n_models + 0.5, "ID")
plt.text(0.05, n_models + 0.5, "Model name")
plt.text(0.4, n_models + 0.5, "Nb members")
plt.text(0.65, n_models + 0.5, "Sea ice model\ncomplexity")
for j_models in range(n_models):
  cplx = cplx_lab[info[j_models][9] - 1]
  plt.text(0, n_models - j_models - 1, str(j_models + 1).zfill(2),               color = colors[j_models])
  plt.text(0.05, n_models - j_models - 1, str(info[j_models][0]),                   color = colors[j_models])
  plt.text(0.4, n_models - j_models - 1, str(len(info[j_models][7])), color = colors[j_models])
  plt.text(0.65, n_models - j_models - 1, cplx,                                     color = colors[j_models])
plt.ylim(1, n_models + 1)
plt.xlim(0, 0.8)
plt.axis("off")
plt.tight_layout()
plt.savefig("./tableED1.png", dpi = 300)
plt.clf()

# Create Fig. S10: CESM ensemble
fig = plt.figure(figsize = (4, 6))
n_members = 35

volume = np.empty(((2100 - 1920 + 1) * 12, n_members))
area   = np.empty(((2100 - 1920 + 1) * 12, n_members))
volume[:] = np.nan
area[:] = np.nan

for j_members in np.arange(1, n_members + 1):
  if j_members == 1:
    ncload("netcdfs/b.e11.B20TRC5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".185001-200512.nc", ["sivol_" + domain, "siarea_" + domain])
    volume1 = eval("sivol_" + domain + "[(1920 - 1850) * 12:]") # subsetting to adjust to other members, that start in 1920
    area1   = eval("siarea_" + domain + "[(1920 - 1850) * 12:]")

    ncload("netcdfs/b.e11.BRCP85C5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".200601-210012.nc", ["sivol_" + domain, "siarea_" + domain])
    volume2 = eval("sivol_" + domain + "[:]")
    area2   = eval("siarea_" + domain + "[:]")
     
  else:
    ncload("netcdfs/b.e11.B20TRC5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".192001-200512.nc", ["sivol_" + domain, "siarea_" + domain])
    volume1 = eval("sivol_" + domain + "[:]")
    area1   = eval("siarea_" + domain + "[:]")
    ncload("netcdfs/b.e11.BRCP85C5CNBDRD.f09_g16." + str(j_members).zfill(3) + ".200601-210012.nc", ["sivol_" + domain, "siarea_" + domain])
    volume2 = eval("sivol_" + domain + "[:]")
    area2   = eval("siarea_" + domain + "[:]")

  volume[:, j_members - 1] = np.concatenate((volume1, volume2), 0)
  area[:, j_members - 1]   = np.concatenate((area1  , area2)  , 0)

yearb = [1920, 1970, 2020]  # All periods that will be tested 
yeare = [1969, 2019, 2069]
colors = [[0.12, 0.56, 1.0], [1.0, 0.56, 0.0], [0.88, 0.08, 0.24]]

for j in range(len(yearb)):
  vm = list()
  nf = list()
  pf = list()

  # 1. Compute the diagnostic for each member
  for j_members in range(n_members):
    nf.append(negative_seaice_feedback(volume[(yearb[j] - 1920) * 12:(yeare[j] - 1920) * 12 + 12, j_members], period = 12)[0])
    pf.append(positive_seaice_feedback(volume[(yearb[j] - 1920) * 12:(yeare[j] - 1920) * 12 + 12, j_members], \
                                       area[  (yearb[j] - 1920) * 12:(yeare[j] - 1920) * 12 + 12, j_members], period = 12)[0])
    vm.append(np.mean(volume[(yearb[j] - 1920) * 12:(yeare[j] - 1920) * 12 + 12, j_members]))

  # 2. Scatter plot
  plt.subplot(2, 1, 1)
  plt.scatter(vm, nf, 8, color = colors[j])
  plt.xlim(-1, 25)
  plt.ylim(-1.0, 0.2)
  plt.grid()
  plt.plot((-1e9, 1e9), (0, 0), 'k-', lw = 2)
  plt.plot((0, 0), (-1e9, 1e9), 'k-', lw = 2)
  plt.ylabel("IFE [m/m]")
  plt.xlabel("Annual mean sea-ice volume [10$^3$ km$^3$]")
  plt.text(np.mean(vm) + 1 + j * (j - 2) * 1.5, np.min(nf) + j * (j - 2) * 0.15, str(yearb[j]) + "-" + str(yeare[j]), color = colors[j])

  plt.subplot(2, 1, 2)
  plt.xlim(-1, 25)
  plt.ylim(-0.2, 1.0)
  plt.grid()
  plt.scatter(vm, pf, 8, color = colors[j])
  plt.plot((-1e9, 1e9), (0, 0), 'k-', lw = 2)
  plt.plot((0, 0), (-1e9, 1e9), 'k-', lw = 2)
  plt.ylabel("OWFE [m$^{-1}$]")
  plt.xlabel("Annual mean sea-ice volume [10$^3$ km$^3$]")
  plt.text(np.max(vm) + 0.5, np.mean(pf), str(yearb[j]) + "-" + str(yeare[j]), color = colors[j])


# Plot CMIP5 envelopes
plt.subplot(2, 1, 1)
xx  = np.arange(1.0, 30.0, 0.05)
uu  = 1.0 / xx
yy1 = (fit[0][0] * uu + fit[0][1])
plt.fill_between(xx, yy1 - 2 * std[0], yy1 + 2 * std[0], linewidth = 0, edgecolor = None, facecolor = [0.8, 0.8, 0.8], zorder = 0)
plt.fill_between(xx, yy1 -     std[0], yy1 +     std[0], linewidth = 0, edgecolor = None, facecolor = [0.7, 0.7, 0.7], zorder = 0)
plt.text(4, -0.15, "CMIP5", color = [0.7, 0.7, 0.7])
plt.text(0.7, 0.1, "a", fontweight = "bold", fontsize = 16)

plt.subplot(2, 1, 2)
yy2 = (fit[1][0] * uu + fit[1][1])
plt.fill_between(xx, yy2 - 2 * std[1], yy2 + 2 * std[1], linewidth = 0, edgecolor = None, facecolor = [0.8, 0.8, 0.8], zorder = 0)
plt.fill_between(xx, yy2 -     std[1], yy2 +     std[1], linewidth = 0, edgecolor = None, facecolor = [0.7, 0.7, 0.7], zorder = 0)
plt.text(8, 0.50, "CMIP5", color = [0.7, 0.7, 0.7])
plt.text(0.7, 0.9, "b", fontweight = "bold", fontsize = 16)


plt.tight_layout()
plt.savefig("./figED5.png", dpi = 300)
