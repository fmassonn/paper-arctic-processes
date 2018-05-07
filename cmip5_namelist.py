# ---------------------------------------
# Namelist defining CMIP5 model structure
# ---------------------------------------
# Notes: CSIRO has sea ice model embedded in atmosphere. Very primitive (probably 0 or 1 layer, no ITD). O'Farrell 1998
#        Russell sea ice model has no ITD
#        SM0L = Semtner 0 layer, no ITD
#        ITDSM0L = Semtner 0 layer, ITD
#        VITD = Virtual ITD
#  
#        Model complexity is going from 1 to 4:
#             Semtner 0-layer   Semtner 3-layer   Winton   Bitz & Lipscomb
#
# no ITD           1                    1            1            1
# virtual ITD      2                    3            3            3
# ITD              3                    4            4            4

#        model name          mod_grid           mask   grid_oce mask_oce yearb   yeare     members                                   #SIM       #Complex

global info

info = [                                                                                    \
        ["ACCESS1-0",        "ACCESS1-0",       100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "CICE4"   , 4], \
        ["ACCESS1-3",        "ACCESS1-3",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "CICE4"   , 4], \
        ["BNU-ESM",          "BNU-ESM",         1.0  ,   "o",     "o",     1850,       2005, [1]                                    ,  "CICE4"   , 4], \
        ["bcc-csm1-1",       "bcc-csm1-1",      100.0,   "o",     "o",     1850,       2012, [1, 2, 3]                              ,  "SIS"     , 4], \
        ["bcc-csm1-1-m",     "bcc-csm1-1",      100.0,   "o",     "o",     1850,       2012, [1, 2, 3]                              ,  "SIS"     , 4], \
        ["CanESM2",          "CanESM2",         100.0,   "a",     "l",     1850,       2005, [1, 2, 3, 4, 5]                        ,  "CanSIM1" , 1], \
        ["CCSM4",            "CCSM4",           100.0,   "o",     "o",     1850,       2005, [1, 2, 3, 4, 5, 6]                     ,  "CICE4"   , 4], \
        ["CESM1-BGC",        "CESM1-BGC",       100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "CICE4"   , 4], \
        ["CESM1-CAM5-1-FV2", "CESM1-BGC",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3, 4]                           ,  "CICE4"   , 4], \
        ["CESM1-CAM5",       "CESM1-BGC",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "CICE4"   , 4], \
        ["CESM1-FASTCHEM",   "CESM1-BGC",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "CICE4"   , 4], \
        ["CESM1-WACCM",      "CESM1-BGC",       100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "CICE4"   , 4], \
        ["CMCC-CMS",         "CMCC-CMS",        100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "LIM2"    , 3], \
        ["CMCC-CM",          "CMCC-CM",         100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "LIM2"    , 3], \
        ["CNRM-CM5",         "CNRM-CM5",        100.0,   "o",     "o",     1850,       2005, [i+1 for i in range(10)]               ,  "GELATO5" , 4], \
        ["CSIRO-Mk3-6-0",    "CSIRO-Mk3-6-0",   100.0,   "a",     "l",     1850,       2005, [i+1 for i in range(10)]               ,  "CSIRO"   , 1], \
        ["EC-EARTH",         "EC-EARTH",        100.0,   "o",     "o",     1850,       2004, [1, 2, 5, 7, 8, 9, 11, 12, 13, 14]     ,  "LIM2"    , 3], \
        ["FGOALS-g2",        "FGOALS-g2",       100.0,   "o",     "o",     1850,       2005, [1, 3, 4, 5]                           ,  "CICE4"   , 4], \
        ["GFDL-CM2p1",       "GFDL-CM3",        1.0,     "o",     "o",     1861,       2005, [i+1 for i in range(9)]                ,  "SIS"     , 4], \
        ["GFDL-CM3",         "GFDL-CM3",        1.0,     "o",     "o",     1860,       2005, [1, 2, 3, 4, 5]                        ,  "SIS"     , 4], \
        ["GFDL-ESM2G",       "GFDL-ESM2G",      100.0,   "o",     "o",     1861,       2005, [1]                                    ,  "SIS"     , 4], \
        ["GFDL-ESM2M",       "GFDL-ESM2M",      1.0,     "o",     "o",     1861,       2005, [1]                                    ,  "SIS"     , 4], \
        ["GISS-E2-H",        "GISS-E2-H",       100.0,   "a",     "l",     1850,       2005, [1, 2, 3, 4, 5, 6]                     ,  "Russell" , 1], \
        ["GISS-E2-R",        "GISS-E2-R",       100.0,   "a",     "l",     1850,       2005, [1, 2, 3, 4, 5, 6]                     ,  "Russell" , 1], \
        ["GISS-E2-H-CC",     "GISS-E2-H",       100.0,   "a",     "l",     1850,       2010, [1]                                    ,  "Russell" , 1], \
        ["GISS-E2-R-CC",     "GISS-E2-R",       100.0,   "a",     "l",     1850,       2010, [1]                                    ,  "Russell" , 1], \
        ["HadCM3",           "HadCM3",          100.0,   "o",     "o",     1860,       2005, [i+1 for i in range(10)]               ,  "SM0L"    , 1], \
        ["HadGEM2-AO",       "HadGEM2-CC",      100.0,   "o",     "o",     1860,       2005, [1]                                    ,  "ITDSM0L" , 3], \
        ["HadGEM2-CC",       "HadGEM2-CC",      100.0,   "o",     "o",     1860,       2004, [1]                                    ,  "ITDSM0L" , 3], \
        ["HadGEM2-ES",       "HadGEM2-CC",      100.0,   "o",     "o",     1860,       2004, [1, 2, 3, 4]                           ,  "ITDSM0L" , 3], \
        ["inmcm4",           "inmcm4",          100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "Unknown" , 3], \
        ["IPSL-CM5A-LR",     "IPSL-CM5A-LR",    100.0,   "o",     "o",     1850,       2005, [1, 2, 3, 4, 5, 6]                     ,  "LIM2"    , 3], \
        ["IPSL-CM5A-MR",     "IPSL-CM5A-MR",    100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "LIM2"    , 3], \
        ["IPSL-CM5B-LR",     "IPSL-CM5A-LR",    100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "LIM2"    , 3], \
        ["MIROC5",           "MIROC5",          1.0,     "o",     "o",     1850,       2012, [1]                                    ,  "ITDSM1L" , 4], \
        ["MIROC-ESM",        "MIROC-ESM",       100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "ITDSM1L" , 4], \
        ["MIROC-ESM-CHEM",   "MIROC-ESM-CHEM",  100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "ITDSM1L" , 4], \
        ["MPI-ESM-LR",       "MPI-ESM-LR",      100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "VITDSM0L", 2], \
        ["MPI-ESM-MR",       "MPI-ESM-MR",      100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "VITDSM0L", 2], \
        ["MPI-ESM-P",        "MPI-ESM-P",       100.0,   "o",     "o",     1850,       2005, [1, 2]                                 ,  "VITDSM0L", 2], \
        ["MRI-CGCM3",        "MRI-CGCM3",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "ITDSM0L?", 4], \
        ["MRI-ESM1",         "MRI-ESM1",        100.0,   "o",     "o",     1851,       2005, [1]                                    ,  "ITDSM0L?", 4], \
        ["NorESM1-M",        "NorESM1-M",       100.0,   "o",     "o",     1850,       2005, [1, 2, 3]                              ,  "CICE4"   , 4], \
        ["NorESM1-ME",       "NorESM1-M",       100.0,   "o",     "o",     1850,       2005, [1]                                    ,  "CICE4"   , 4], \
       ]
