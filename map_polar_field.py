#!/usr/bin/python
# Author: Francois Massonnet
# Date:   Jan 2017
# Contents:
# 1) A function that plots fields in polar stereo projection
#    The input to this function is a numpy array, not a 
#    NetCDF file. This function can be called from within a script,
#    not from the shell command
#    
# 2) A script (below) that takes NetCDF files, and then calls the function 1)
#    This script makes the link with the shell

# ------------------------------
# PART 1) The plotting function |
# ------------------------------
def map_polar_field(field                      , \
                    lon         = None         , \
                    lat         = None         , \
                    hemisphere  = "n"          , \
                    lb          = None         , \
                    ub          = None         , \
                    vname       = "unknown var", \
                    colmap      = "RdYlBu_r"   , \
                    nsteps      = 20           , \
                    extmethod   = "both"       , \
                    units       = "units"      , \
                    plain_continents = True    , \
                    title       = "Title"      , \
                    figdir      = ""           , \
                    filename    = "map_field"  , \
                    imageformat = "png"        , \
                    dpi         = 90           ):
  """ Purpose:      Plots a 2-D horizontal field in polar projection

       - field      is a 2-D numpy array. Since the grid can be irregular,
                    nothing is assumed about the order of dimensions

       - lon        is a 1-D or a 2-D numpy array, representing longitudes. 
                    If lon is 1-D, it is expected to describe a regular grid.
                    If None, a global regular grid is assumed and lon is built
                    to describe this grid.

       - lat        is a 1-D or a 2-D numpy array, representing latitudes.
                    If lat is 1-D, it is expected to describe a regular grid.
                    if None, a global regular grid is assumed and lat is built
                    to dsecribe this grid.

       - hemisphere is the hemisphere to plot.

       - lb, ub     are lower and upper bounds for the colorbar. If None, they
                    are respectively set to be the minimum and maximum of the 
                    field over the plotted domain

       - vname      is the name of the variable plotted

       - colmap     is the colormap to be used

       - nsteps     is an indicative number of steps to put on the colorbar.
                    This number will slightly change in the final plot in order
                    to have readable steps.

       - extmethod  is the type of color bar (<-->; |-->; <--|; |--|)

       - units      is the string describing the units of the variable plotted
  
       - plain_continents is a boolean to be set to True if continents are to
                    be filled

       - figdir     is the output directory where the figure is to be printed

       - title      is the title to display 

       - filename   is the name of the figure

       - imageformat is the format of the figure (eps, png, pdf)

       - dpi        is the resolution of the figure in dots per inch.


      Author:   Francois Massonnet
      Date  :   February 1st 2017

  """

  # Import necessary libraries and modules
  # --------------------------------------
  import numpy             as np
  import matplotlib; matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  import sys
  from mpl_toolkits.basemap import Basemap, addcyclic

  # ---------------------------------------------
  # 1. Check that the data is coming as expected |
  # ---------------------------------------------

  # Check that field has exactly two dimensions
  if len(field.shape) != 2:
    sys.exit("(map_polar_field): field has not exactly two dimensions!")

  # Check that longitude and latitude don't have different dimensions
  if (lon is not None and lat is not None) and \
     (len(lon.shape) != len(lat.shape)):
    sys.exit("(map_polar_field): longitudes and latitudes dimension differ")

  # -----------------------------
  # 2. Retrieve grid information |
  # -----------------------------

  # Call them arbitrarily nx and ny (1st and 2nd)
  nx = field.shape[0]
  ny = field.shape[1]

  # If lon and/or lat aren't defined, then we assume a regular global grid 
  if lat is None:
    lat = np.linspace(-90, 90, nx)
  if lon is None:
    lon = np.linspace(0, 360, ny)

  # Add cyclicity in case of 1-D longitude (means regular grid)
  if len(lon.shape) == 1:
    field, lon = addcyclic(field, lon)

  # If lon and lat are 1-D, mesh the domain 
  # Note that if lat is 1-D then lon is 1-D too automatically (see above tests)
  if len(lat.shape) == 1:
    lon, lat = np.meshgrid(lon, lat)

  # -----------------------
  # 3. Create the colorbar |
  # -----------------------

  # If lower and upper bounds aren't defined, determine them from the field
  if lb is None:
    lb = np.max([1e9, np.min(field)])
  if ub is None:
    ub = np.max(field)
  # Create a colormap so that it has readable numbers. This implies having
  # a number of steps close to the one provided, and intervals that are not
  # crazy numbers. We want to avoid this type of situation:
  # [0.0, 0.111111, 0.222222, 0.333333, ...]

  step = 1.0
  notconverged = True
  ntry = 0
  while notconverged and ntry <= 100:
    ntry += 1
    if len(np.arange(lb, ub + step, step)) <= nsteps * 0.5:
      step = step / np.max((2.0, (-1) ** ntry * 5.0))
    elif len(np.arange(lb, ub + step, step)) >= nsteps * 1.5:
      step = step * np.max((2.0, (-1) ** ntry * 5.0))
    else:
      notconverged = False

  if notconverged:
    print("(map_polar_field). Unable to provide a colorbar")
    print(" Default colorbar will be used")
    clevs = np.arange(lb, ub, nsteps)
  else: 
    clevs = np.arange(lb, ub + step, step)

  # Load the colormap
  cmap = eval("plt.cm." + colmap)

  # Colors for values outside the colorbar
  # Get first and last colors of the colormap

  first_col = cmap(0)[0:3]
  last_col  = cmap(255)[0:3]
  # Create "over" and "under" colors.
  # Take respectively the latest and first color and
  # darken them by 50%
  col_under = [i / 2.0 for i in first_col]
  col_over  = [i / 2.0 for i in last_col ]
 
  # -----------------------
  # 4. Create a projection |
  # -----------------------

  # Determine if we are in the northern hemisphere or in the Southern hemisphere.
  if hemisphere == "n":
    boundlat = 60.
    l0 = 0.
  elif hemisphere == "s":
    boundlat = -50.
    l0 = 180.
  else:
    sys.exit("(map_polar_field) Hemisphere unkown")

  # Projection name
  projname = hemisphere + "plaea"
  
  # Initialize projection (largely inspired from
  # http://matplotlib.org/basemap/users/examples.html

  map = Basemap(projection = projname, boundinglat = boundlat, \
                lon_0 = l0, resolution = 'l')

  # Map the lon, lat to the projection
  x, y = map(lon, lat)

  # -------------------
  # 5. Plot the figure |
  # -------------------

  # Open figure
  plt.figure("fig", figsize = (6, 6))

  # Draw coastlines, country boundaries, fill continents, meridians & parallels
  if plain_continents:
    map.drawcoastlines(linewidth = 0.25)
    map.fillcontinents(color = 'grey', lake_color = 'w')
  else:
    map.drawcoastlines(linewidth = 1.0)

  map.drawmeridians(np.arange(0, 360, 30))
  map.drawparallels(np.arange(-90, 90, 10))
    

  # A feature of the NEMO Grid is that Basemap produces super-strange results
  # if the fields x, y, field are plotted
  # Plotting them up to the last line fixes it.
  if hemisphere == "s":
    x = x[0:-1, 0:-1]
    y = y[0:-1, 0:-1]
    field =  field[0:-1, 0:-1]

  # Create a contourf object called "cs" 
  cs = map.contourf(x, y, field, clevs, cmap = cmap, \
                    latlon = False, extend = extmethod)

  cs.cmap.set_under(col_under)
  cs.cmap.set_over(col_over)

  cbar = map.colorbar(cs, location = 'bottom', pad = "5%")

  cbar.set_label(units)

  xt = cbar.ax.get_xticklabels()
  cbar.ax.set_xticklabels(xt, rotation = 0)

  plt.title(title)

  # Save figure
  plt.savefig(figdir + filename + "." + imageformat, bbox_inches = "tight", dpi = dpi)
  print('Figure '+ figdir + filename + "." + imageformat + " printed")
  plt.close("fig")





# -------------------------------------------------------------
# PART 2) The script to communicate with outside world (shell) |
# -------------------------------------------------------------
if __name__ == '__main__':
  
  from netCDF4 import Dataset
  import numpy as np
  import sys
  import os
  
  def help():
    print("")
    print("  Usage: map_polar_field.py NETCDF-file VARIABLE-name")
    print("      PURPOSE : ")
    print("        Makes the maps (polar stereographic projections) of the")
    print("        variable passed as argument")
    print("      ARGUMENTS : ")
    print("        NETCDF-file : NetCDF file in which the field to be plotted")
    print("                      can be found")
    print("        VARIABLE-name :")
    print("                      Name of the variable to be plotted (must be in")
    print("                      the NetCDF file")
    print("        LOWER-bound : Lower bound to plot")
    print("        UPPER-bound : Upper bound to plot")
    print("      HELPDESK : ")
    print("        francois.massonnet@bsc.es")
    print("")
  
  
  # --------------------------------------------------------------
  # 1. Read the number of arguments, and the arguments themselves |
  # --------------------------------------------------------------
  
  # If not five arguments are passed (script name + file name + variable + lower bound + upper bound), abort
  if len(sys.argv) <= 1:
    help()
    sys.exit("No argument given, two expected.")
  elif len(sys.argv) == 2:
    filein = sys.argv[1]
    if not os.path.isfile(filein):
      sys.exit("File " + filein + " does not exist")
    else:
      f = Dataset(filein, mode = "r")
      print("One argument given, four expected.")
      print("Second argument must be one of the following variables:")
      a = [v.encode('ascii') for v in f.variables]
      print(a)
      f.close()
      sys.exit()
  elif len(sys.argv) != 5:
    sys.exit("Not four arguments given")
  
  
  # Retrieve arguments
  filein = sys.argv[1]
  varin  = sys.argv[2]
  lb     = np.float32(sys.argv[3])
  ub     = np.float32(sys.argv[4])

  # ---------------------------------------------
  # 2. Retrieve information from the NetCDF file |
  # ---------------------------------------------
  
  # Open NetCDF and retrieve variable and its units
  f = Dataset(filein, mode = "r")
  
  try:
    field = f.variables[varin][:]
  except KeyError:
    f.close()
    sys.exit("Variable " + varin + " not found in the NetCDF file provided")
  
  try:
    units = f.variables[varin].units
  except AttributeError:
    units = "units not available"
  # Determine the number of dimensions of the field.
  # If there are two dimensions, this is is simple (they are assumed
  # to be the spatial dimensions)
  nt = 1
  nz = 1
  isTime = False
  isDepth = False
  
  if len(field.shape) < 2:
    sys.exit("The input variable has less than two dimensions")
  elif len(field.shape) == 2:
    # If there are two dimensions, name them arbitrarily nx and ny
    # (remember that they are not necessarily lon and lat, since
    # grids may be irregular.
    nx = field.shape[0]
    ny = field.shape[1]
  else:
    # If more than 2 dimensions, let's look for these other dimensions.
    # Identify whether there are time or other dimensions
    # than spatial
    for v in f.variables[varin].dimensions:
      if v == "x" or v == "y" or v == "i" or v == "j" or v == "lat" or v == "lon" \
         or v == "longitude" or v == "latitude":
        pass # skip, these are known!
      elif v == "time" or v == "time_counter" or v == "t":
        isTime = True
        nt = f.dimensions[v].size
      elif v == "z" or v == "depth" or v == "deptht" or v == "level" or v == "lev":
        isDepth = True
        nz = f.dimensions[v].size
      else:
        sys.exit("Dimension " + v + " not referenced")
  
  # Retrieve information about longitude and latitude
  lonnotfound = True
  lonlist = ["lon", "long", "longitude", "longitudes", "glamt", \
             "nav_lon", "nav_lon_grid_T"]
  j = 0
  while lonnotfound and j < len(lonlist):
    try:
      lon = f.variables[lonlist[j]][:]
      lonnotfound = False
      break
    except KeyError:
      j += 1
  
  # If still nothing is found about longitude, look into the dimensions
  if lonnotfound:
    for d in f.variables[varin].dimensions:
      if d == "longitude" or d == "lon":
        nlon = f.dimensions[d].size
        lon = np.linspace(0, 360, nlon)
  
  # Repeat for lat
  latnotfound = True
  latlist = ["lat", "latitude", "latitudes", "gphit", \
             "nav_lat", "nav_lat_grid_T"]
  j = 0
  while latnotfound and j < len(latlist):
    try:
      lat = f.variables[latlist[j]][:]
      latnotfound = False
      break
    except KeyError:
      j += 1
  
  if latnotfound:
    for d in f.variables[varin].dimensions:
      if d == "latitude" or d == "lat":
        nlat = f.dimensions[d].size
        lat = np.linspace(-90, 90, nlat)
  
  
  # ----------------------------------
  # 3. Determine graphical properties |
  # ----------------------------------
  # We say that a symmetric colorbar should be used
  # if the average of the field is less than 1% of the
  # max (arbitrary)
  if -lb == ub:
    colmap = "RdBu_r"
    # Determine lower and upper bounds automatically
    extmethod = "both"
  else:    
    # Else, we browse through a catalogue of possible variables
    if varin == "siconc" or varin == "sic" or varin == "ice_conc":
      colmap = "Blues_r"
      extmethod = "neither"
    elif varin == "sivolu" or varin == "sithick" \
      or varin == "sithic" or varin == "sivol":
      colmap = "terrain"
      extmethod = "max"
    elif varin == "sst" or varin == "tos":
      colmap = "CMRmap"
      extmethod = "both"
    elif varin == "sss" or varin == "sos" or varin == "so":
      colmap = "winter"
      extmethod = "both"
    elif varin == "tas" or varin == "t2m":
      colmap = "plasma"
      extmethod = "both"
    elif varin == "correlation":
      colmap = "RdBu_r"
      extmethod = "neither"
    else:
      colmap = "cubehelix_r"
      extmethod = "both"

      # Determine colorbar in a user-friendly way. 
      # A first guess of the lower and upper bounds is determined
      # from the 5th and 95th percentiles.
      lbt = lb
      ubt = ub
      if lbt == 0.0:
        extmethod = "max"
        lb = 0.0
      else:
        # Now it can happen that the value of lbt is not convenient, 
        # say 1.3467. We would like to round that to readable values.
        # This is what this hack does (also works for 0.000XYZ numbers)
        si = np.sign(lbt) # The sign of the lower bound
        ex = np.floor(np.log10(np.abs(lbt))) # The order (power of ten)
                                             # ex. 0.001 --> -3

        # The rounded leading digit. Ex 0.000432 --> 0.4; 12.431 --> 12.
        pa = np.floor(np.abs(lbt) / (10.0 ** np.floor(np.log10(np.abs(lbt)))))

        # Re-definining the lower bound, now rounded         
        lb = si * pa * 10.0 ** ex
      if ubt == 0.0:
        extmethod = "min"
        ub = 0.0
      else:
        # See above, this is the same procedure
        si = np.sign(ubt)
        ex = np.floor(np.log10(np.abs(ubt)))
        pa = np.ceil(np.abs(ubt) / (10.0 ** np.floor(np.log10(np.abs(ubt)))))
        
        ub = si * pa * 10.0 ** ex

  # Determine hemisphere(s) to plot
  if np.max(lat) < 0:
    hemi = ("s")
  elif np.min(lat) > 0:
    hemi = ("n")
  else:
    hemi = ("n", "s")
  
  # Define a figure name
  fname = filein.split("/")[-1] 
  # -------------------------------------
  # 4. Call the map_polar_field function |
  # -------------------------------------
  for jt in range(nt): 
    for jz in range(nz):
      for hemisphere in hemi:
        fname = filein.split("/")[-1] + "_" + varin + "_" + hemisphere + "_t" + \
                str(jt + 1).zfill(3) + "_z" + str(jz + 1).zfill(3) 
        title = filein.split("/")[-1] + "\n" + varin + ": t = " + \
                str(jt + 1).zfill(3) + ", z = " + str(jz + 1).zfill(3)
        if isTime and not isDepth:
          field_to_plot = np.squeeze(field[jt, :, :])
        elif not isTime and isDepth:
          field_to_plot = np.squeeze(field[jz, :, :])
        elif isTime and isDepth:
          field_to_plot = np.squeeze(field[jt, jz, :, :])
        else:
          field_to_plot = field
  
        map_polar_field(field_to_plot               , \
                         lon         = lon          , \
                         lat         = lat          , \
                         hemisphere  = hemisphere   , \
                         lb          = lb           , \
                         ub          = ub           , \
                         vname       = varin        , \
                         colmap      = colmap       , \
                         extmethod   = extmethod    , \
                         units       = units        , \
                         title       = title        , \
                         filename    = fname        , \
                                                    )
