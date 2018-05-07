#!/usr/bin/python
#
# Creation of Custom colormaps for Python
from matplotlib.colors import LinearSegmentedColormap


cdict  = {'red':   (
                    (0.00, 0.00, 0.00)    ,
                    (0.60, 0.86, 0.86)    ,
                    (1.00, 0.30, 1.00))   ,

          'green': ((0.00, 0.22, 0.22)    ,
                    (0.60, 0.86, 0.86)    , 
                    (1.00, 0.00, 0.58))   ,

          'blue':  ((0.00, 0.35, 0.35)    ,
                    (0.60, 1.00, 1.00)    ,
                    (1.00, 0.64, 0.64)    , 
                   )
        }


global mycmap
mycmap = LinearSegmentedColormap('mycmap', cdict)

#######################
# Sea ice concentration
#######################

cdict  = {'red':   (
                    (0.00, 0.00, 0.00)    ,
                    (0.60, 0.50, 0.50)    ,
                    (0.80, 0.10, 0.10)    ,
                    (1.00, 1.00, 1.00))   ,

          'green': ((0.00, 0.00, 0.00)    ,
                    (0.60, 0.50, 0.50)    ,
                    (0.80, 1.00, 1.00)    ,
                    (1.00, 1.00, 1.00))   ,

          'blue':  ((0.00, 0.40, 0.40)    ,
                    (0.60, 0.95, 0.95)    ,
                    (0.80, 1.00, 1.00)    ,
                    (1.00, 1.00, 1.00)    ,
                   )
        }


global mycmap_sic
mycmap_sic = LinearSegmentedColormap('mycmap_sic', cdict)

####################

