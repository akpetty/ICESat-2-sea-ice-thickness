# Sea ice thickness with ICESat-2
### Project team 
Alek Petty, Nathan Kurtz, Ron Kwok, Thorsten Markus, Tom Neumann   

Get in touch if you have any questions:
[alek.a.petty@nasa.gov](alek.a.petty@nasa.gov) / [alekpetty.com](http://www.alekpetty.com) / [@alekpetty](https://twitter.com/alekpetty)   

## Introduction

Code repository for producing sea ice thickness estimates from ICESat-2 freeboard data. The plan is for this repository to include all the code and methodology needed for converting sea ice freeboard to sea ice thickness from data released after the upcoming launch of ICESat-2. This will include testing with ICESat-1 and Operation IceBridge freeboard data. Our primary focus will be on the Arctic, but we do plan on including Southern Ocean thickness estimates at some point. Jupyter Notebooks included in this repo will provide an interactive view of the freeboard-thickness methodology.

![Sea ice thickness from satelite laser altimetry](SeaIceSchematic.png?raw=true "Sea ice thickness from satelite laser altimetry")


The main inputs needed for the thickness calculation are sea ice freeboard and snow depth/density. We also need to make a few other assumptions (e.g. ice density) which are discussed in the scripts/notebooks.

Included in this repository are some test along-track freeboard data files from:
 - ICESat-1 (??): ```/TestData/IS1/```
 - Operation IceBridge (Apr 25th, 2013): ```/TestData/OIB/```
 - ICESat-2 data (one file of ATL10 along-track shot freeboards, Oct 14th 2018): ```/TestData/IS2/```
and gridded snow depths/densities from the NESOSIM model:
- NESOSIM (Aug 15th to May 1st 2005-2006, 2013-2014, Oct 2018): ```/TestData/NESOSIM/```

The ICESat-2 ATL10 along-track freeboard data is being made available on the SCF ('/cooler/sea_ice_pge/sea_ice_product/') so the plan is to run this code on that server and link directly to those data files as they are produced. 
 
More information about the NESOSIM snow model can be found here: https://github.com/akpetty/NESOSIM. I also have a working development branch (NESOSIMdev) that I can add you to as desired. I have included a test file run through to the end of October 2018 to use with the October ATL10s.

The ```/Scripts/``` folder contains some Jupyter Notebooks and Python scripts describing the snow distribution and thickness conversion methodology.


I will include more information on getting started with Python etc. The plan is for the code to be produced in Python 3.6. Packages will be maintained with Conda. You can get more information about this in my sea ice prediction repo: https://github.com/akpetty/SeaIcePrediction. Briefly, you should be able to just import the conda environment included in this repo as:
```
conda env create -n py36env -f environment.yml
```
then simply activate this environment to work in (to avoid using you own Python environment which wil likely have dependency issues)
```
source activate py36env
```

## Project plan

![Sea ice thickness flowchart](IS2flowchart.pdf?raw=true "Sea ice thickness processing flowchart")


## Data description

Add here

## Script/notebook descriptions

Add here

## References

Add references
