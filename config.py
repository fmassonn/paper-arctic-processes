# Set of modules to be systematically imported
# for all Python scripts

# ---------------------------------------
# Headers, modules imported, etc.
# ---------------------------------------
# Regular Python modules
import matplotlib
matplotlib.use('Agg') # avoids tKinter error
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.dates as mdates

from mpl_toolkits.basemap import Basemap, addcyclic

import csv
import numpy as np
import sys
import datetime

from netCDF4 import Dataset
from matplotlib.gridspec import GridSpec

import os

# Change font globally
# --------------------
font_dirs  = [u'./fonts/', ]
font_files = font_manager.findSystemFonts(fontpaths = font_dirs)
font_list  = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)
matplotlib.rcParams['font.family'] = 'Acephimere'

# In-home functions
# -----------------
exec(open("ncload.py").read())     # To read NetCDF quickly
from   seaice_commondiags import * # To compute basic sea ice diagnostics
from   map_polar_field    import * # To map fields with polar projection
from   ts_analyses        import * # To conduct time series analyses (detrending, e-folding time scale, ...)

# In-home colormaps
# -----------------
exec(open("custom_colormaps.py").read())

plt.close("all")                   # Close all figures

