#!/usr/bin/python
#
# Function to load NetCDF file in a convenient way
# F. Massonnet
#

def ncload(filein, variables = None):
  """ FUNCTION ncload
  
      USAGE:
        ncload(filein, variables = None, times = None):
      
      PURPOSE:
        Reads a NetCDF file and loads the data it contains 
        (or a part of it) in the workspace

      ARGUMENTS:
        (mandatory) filein    [string] 
          Path to the NetCDF file
        
        (optional)  variables [list of strings]
          The variables to be loaded
          If variables is unspecified or set to None, all variables
          found will be loaded

      EXAMPLE:
        Suppose a NetCDF file named "file.nc" located in the directory
        /path/to/file  has three variables "lon", "lat", "temperature"

        ncload("/path/to/file.nc")
          will load "lon", "lat", "temperature"
         
        ncload("/path/to/file.nc", ["lat", "temperature"])
          will only load lat and temperature

      REMARKS:
        There seems to be problems when two NetCDF files with the same
        variable names are read successively.

      CONTACT:
        Francois Massonnet, francois.massonnet@uclouvain.be
  """

  from netCDF4 import Dataset
  import numpy as np
  import sys
  
  try:
    f = Dataset(filein, mode = "r")
  except IOError:
    sys.exit("(ncload) ERROR: File not found:\n >>> " + (filein) + " <<<")

  # List of variables
  if variables is None:
    variables = [v.encode('ascii') for v in f.variables]
  else:
    if type(variables) is not list:
      if type(variables) is str: # convert to list if only one variable given
        variables = [variables]
      else:
        sys.exit("(ncload) ERROR: input variables not provided as a list")
    
  l = list()

  for v in variables:
    print("(ncload) -- loading " + str(v).ljust(20) + " in environment")
    exec("try:\n"                        +\
         "  global " + str(v) + "\n"     +\
         "  " + str(v) + " = np.squeeze(f.variables[\"" + v + "\"][:])\n" \
         "except KeyError:\n"            +\
         "  sys.exit(\"(ncload) ERROR: variable \" + str(v) + \" of file \" + filein + \" does not exist in NetCDF file\")"                              )
    
  f.close()
