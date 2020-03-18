# Sea ice thickness with ICESat-2
Contact: Alek Petty / alek.a.petty@nasa.gov / www.alekpetty.com

Code repository for producing sea ice thickness estimates from ICESat-2 freeboard data. The model code is written in the open source programming language Python (Python Software Foundation, https://www.python.org/). The main inputs needed for the thickness calculation are ATL10 sea ice freeboards from ICESat-2 and snow depth/density data (e.g. from the NESOSIM model). We also need to make a few other assumptions (e.g. ice density) which are discussed in the code.

![Sea ice thickness from satelite laser altimetry](SeaIceSchematic.png?raw=true "Sea ice thickness from satelite laser altimetry"{:height="70%" width="70%"})

Versions:

v1.0: This initial repository was used to processes ICESat-2 (and also ICESat) thickness estimates for the first winter season of data collection (October 14th 2018 to April 30th 2019).

### Getting Started

If you are familiar with using conda to manage your Python environment you should be able to just import the conda environment included in this repo as:
```
conda env create -n py36envX -f environment.yml
```
then simply activate this environment to work in (to avoid using you own Python environment which wil likely have dependency issues)
```
source activate py36envX
```
A requirements.txt file is also included for installing with pip, although this hasn't been tested extensively.


### Code layout

The ```/Code/``` folder contains the primary scripts to process ATL10 freeboards, and also ICESat freeboards, into sea ice thickness. The ```batch_process_icesat2.py``` and ```batch_process_icesat.py``` file links to modules contained in ```common_functions.py```

Also included are plotting scripts:
```/Code/Plotting/``` 
which include code to generate some of the primary plots shown in Petty et al., (2020). 
 
![Sea ice thickness flowchart](IS2flowchart.png?raw=true "Sea ice thickness processing flowchart"{:height="70%" width="70%"})

## Input data

#### ATL10 freeboards

We use the ICESat-2 ATL10 sea ice freeboard product (designated Release 002) which is disseminated through the NSIDC (Kwok et al., 2019b, https://nsidc.org/data/atl10). The six beams are comprised of three beam pairs, with each beam pair containing a strong and weak beam which are separated by 90 m across-track and 2.5 km along-track, with each beam pair then separated by around 3.3 km across-track. Individual segment heights are produced from each beam using 150-photon aggregates, in an effort to produce heights with a precision of 2 cm or less over flat surfaces, as described in the ATL07 sea ice/sea surface height product description and ATBD (Kwok et al., 2019a; Kwok et al., 2019c). This results in segment lengths of around 10 m to 200 m for the strong beam (mean of around 15 m), and around 40 m to 800 m (mean of around 60 m) for the weak beam (Kwok et al., 2019c). We primarily utilize strong beam 1 for our processing but have carried out beam consistency checks to show this has only a small impact on our basin-scale distributions. 

#### NESOSIM snow depth and density

We primarily make use of snow depth and density data from the NASA Eulerian Snow On Sea Ice Model (NESOSIM) v1.0; a new open-source snow budget model that is currently configured to simulate snow on sea ice across the Arctic Ocean through the accumulation season (Petty et al., 2018). More information about the NESOSIM snow model can be found here: https://github.com/akpetty/NESOSIM. For this time period we run NESOSIM using ERA-Interim (ERA-I) snowfall and winds (Dee et al., 2011), NASA Climate Data Record (CDR) sea ice concentrations (Meier et al., 2017) and the European Organization for the Exploitation of Meteorological Satellites (EUMETSAT) Ocean and Sea Ice Satellite Application Facility (OSI SAF, www.osi-saf.org) ice drifts (Lavergne et al., 2010).  

#### OSI SAF ice type

We use the EUMETSAT OSI SAF sea ice type product (Breivik et al., 2012) to delineate our results between FYI and MYI. This is used mainly to derive modified Warren climatology snow depth information, but is also stored in the datafiles to enable post-processing. The data can be accessed at (http://www.osi-saf.org/?q=content/global-sea-ice-type-c)

#### Ancillary data

We also utilize the NSIDC regional mask of the Arctic Ocean and its peripheral seas to delineate the results by Arctic region (region_n.msk in /AncData/). 

## Derived data

Data are in the process of being provided on the NSIDC. For now just the gridded data are available on the ICESat-2 project page: PROVIDE LINK.

#### Along track data


#### Gridded data

A single NetCDF file is provided for each month from November 2018 to April 2019. The plan is to provide this routinely for all future winter seasons (just for the months of October through to April when we have higher confidence of the snow loading and reduced likelihood of surface melt complexity).

Each NetCDF file includes the following (monthly mean) variables: 

 - Sea ice freeboard (meters)
 - NESOSIM snow depth (meters)
 - Sea ice thickness (meters)
 - Sea ice thickness uncertainty (meters)
 - NESOSIM snow density (kilograms per meter cubed)
 - OSI SAF sea ice type (meters)
 - Effective day of the month (day)

![April 2019 gridded sea ice thickness dataset](april2019_gridded_demo.png?raw=true "April 2019 gridded sea ice thickness dataset"{:height="70%" width="70%"})


## Reference

Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann, Winter Arctic sea ice thickness from ICESat-2 freeboards (JGR Oceans)


