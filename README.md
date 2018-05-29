# README

This repository contains the scripts and data necessary to replicate bit wise the results of the following paper:

Massonnet, F. et al. Arctic sea-ice change tied to its mean state through thermodynamic processes

accepted in Nature Climate Change, doi:10.1038/s41558-018-0204-z

Any question about the use of this repository should be directed to the corresponding author, François Massonnet (francois.massonnet@uclouvain.be)

## Getting the data

In order to run the scripts of this project, data must be first retrieved. 

Data can be retrieved from the PANGAEA Data Publisher under the following reference:

https://doi.pangaea.de/10.1594/PANGAEA.889757

The data is available as a ZIP file named `SeaIce_model_data.zip`. In the main directory, create a directory named `netcdfs`, move the zip file there, and unzip it. After this operation, the structure of the folder therefore looks like:

```
paper-arctic-processes/
├── *.py        # All Python scripts
├── csv/
├── fonts/
├── netcdfs/
│   ├── *.nc    # NetCDF files downloaded from the archive
├── README.md   # The current file
```

## Running the scripts

Scripts are self-documented. They are written in Python. Python2.7 will certainly work, the scripts have not been tested with Python3 but could work as well.

