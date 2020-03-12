# Sea ice thickness with ICESat-2
### Project team 
Alek Petty, Nathan Kurtz, Ron Kwok, Thorsten Markus, Tom Neumann   

Get in touch if you have any questions:
[alek.a.petty@nasa.gov](alek.a.petty@nasa.gov) / [alekpetty.com](http://www.alekpetty.com) / [@alekpetty](https://twitter.com/alekpetty)   

## Introduction

Code repository for producing sea ice thickness estimates from ICESat-2 freeboard data. The model code is written in Python, an open source programming language (Python Software Foundation, https://www.python.org/). The main inputs needed for the thickness calculation are ATL10 sea ice freeboard from ICESat-2 and snow depth/density. We also need to make a few other assumptions (e.g. ice density) which are discussed in the code.

![Sea ice thickness from satelite laser altimetry](SeaIceSchematic.png?raw=true "Sea ice thickness from satelite laser altimetry")

Versions:

v1.0: This initial repo was used to processes ICESat-2 (and also ICESat) thickness estimates for the first winter season of data collection (October 14th 2018 to April 30th 2019).

### Getting Started
If you are familiar with using conda to manage your Python environment you should be able to just import the conda environment included in this repo as:
```
conda env create -n py36envX -f environment.yml
```
then simply activate this environment to work in (to avoid using you own Python environment which wil likely have dependency issues)
```
source activate py36envX
```
A requirements.txt file is also included for installing with pip, although this hasn't been tested.



## Code layout


The ```/Code/``` folder contains the primary scripts to process ATL10 freeboards, and also ICESat freeboards into sea ice thickness. The ```batch_process_icesat2.py``` file links to modules contained in ```common_functions.py```

Also included are plotting scripts:
```/Code/Plotting/``` 
which includes scripts to generate the plots shown in Petty et al., (2020). 

We also include gridding scripts (for both the ICESat and ICESat-2 data):
```/Code/Gridding/``` 

![Sea ice thickness flowchart](IS2flowchart.png?raw=true "Sea ice thickness processing flowchart")

## Input data

#### ATL10 freeboards

#### NESOSIM snow depth and density

More information about the NESOSIM snow model can be found here: https://github.com/akpetty/NESOSIM. 

#### OSI SAF ice type

This is used mainly to derive modified Warren climatology snow depth information, but is also stored in the datafiles to enable post-processing.

#### Ancillary data

## Derived data

Data are in the process of being provided on the NSIDC. For now just the gridded data are available on the ICESat-2 project page: PROVIDE LINK.

#### Along track data


#### Gridded data

A single NetCDF file is provided for each month from November 2018 to APril 2019. The plan is to provide this routinely for all future winter seasons (October to April).
Each NetCDF file includes the following variables: 

 - Effective day of the month (day). 

 PROVIDE DEMO MONTH FILE HERE.

## Reference

Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann, Winter Arctic sea ice thickness from ICESat-2 freeboards (JGR Oceans)


