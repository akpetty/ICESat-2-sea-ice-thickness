""" BatchProcessIS2shot.py
	
	Processing sea ice thickness with ICESat-2 freeboards
	Initial code written by Alek Petty (01/06/2019)
	
	Input:
		ICESat-2 freeboards, NESOSIM snow depths/densities, OSISAF ice type

	Output:
		Ice thickness using various snow depth assumptions.

	Python dependencies:
		See below for the relevant module imports. Of note:
		xarray/pandas
		netCDF4
		matplotlib
		basemap

	GSFC: W99mod5 snow depth/density, 915 kg/m3 ice density
	CPOM: regional mean W99mod5, 882, 917 kg/m3 ice density
	JPL: W99mod5 (bus using ice type fraction), 
	AWI: W99mod5, 882, 917 kg/m3 ice density
	Petty: NESOSIM distributed, 899 MYI, 916 

		More information on installation is given in the README file.

	Update history:
		01/06/2019: Version 1.
    
"""

import matplotlib, sys
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
from glob import glob
import netCDF4 as nc4
from scipy.interpolate import griddata
from netCDF4 import Dataset

import time
import sys
import common_functions as cF

from itertools import repeat, product
import concurrent.futures

#import multiprocessing
#from numba import jit, prange


def main(freeboardFile):
	""" Main ICESat-1 processing module 
	
	Convert the ATL10 shot freeboard data to ice thickness using various snow loading input estimates
	
	# add distribute modified Warren snow depths.
	# add Ron snow depths. 
	
	Args:
		fileT (str): the ATL10 file path

	"""

	
	outStr=freeboardFile.split("/")[-1][:-3]

	regionflags=1
	icetype=0
	warrensnow=1
	modwarrensnow5=0
	nesosim=1
	nesosimdisttributed=1
	savedata=0
	saveNetCDFX=1
	
	
	# Map projection
	mapProj = Basemap(projection='npstere',boundinglat=55,lon_0=0, resolution='l' , round=False)

	# Get freeboard files
	start = time.time()
	print ('Grabbing freeboard data...')
	dF = cF.getIS1FreeboardData(freeboardFile, mapProj)
	print ('Got freeboard data in '+str(np.round((time.time()-start), 2))+' seconds')

	if isinstance(dF, pd.DataFrame):
		print('Good file') 
	else:
		# Bad file so skip
		print('Bad file')
		return

	if (dF.shape[0]<10):
		# If file had no good freeboard data break
		print('Not enough data in file', dF.count)
		return
	print ('Got good ATL10 freeboard data in '+str(np.round((time.time()-start), 2))+' seconds')
	

	# ----- Region flags -------
	if (regionflags==1):
		print ('Assign NSIDC region mask...')
		dF=cF.assignRegionMask(dF, mapProj, ancDataPath=ancDataPath)
		#print(dF.head(3)) 
		print ('Complete')

	
	# ----- Get colocated ice type data -------
	if (icetype==1):
		start = time.time()
		print ('Processing ice type data...')
		dF = cF.getIceType(dF, iceTypePath, mapProj, res=4, returnRaw=0)
		print ('Processed ice type in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
	
	#------- Get Warren snow depths -----------
	if (warrensnow==1):
		start = time.time()
		print ('Processing Warren snow depths...')
		dF = cF.getWarrenData(dF, outSnowVar='snow_depth_W99', outDensityVar='snow_density_W99')
		print ('Processed Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
	
	#------- Get modified (50%) Warren snow depths -----------
	if (modwarrensnow5==1):
		start = time.time()
		print ('Processing Warren mod5 snow depths...')
		dF = cF.getWarrenData(dF, outSnowVar='snow_depth_W99mod5', outDensityVar='None', modFactor=0.5)
		print ('Processed Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)

	#------- Get NESOSIM snow depths -----------
	if (nesosim==1):
		# ----- Get dates and coincident NESOSIM file -------
		print ('Grabbing NESOSIM file and dates...')
		dateStr= cF.getDate(dF['year'].iloc[0], dF['month'].iloc[0], dF['day'].iloc[0])
		print('NESOSIM date:', dateStr)
		fileNESOSIM, dateStr = cF.getNesosimDates(dF, snowPath)
		print ('Complete')
		#print ('NESOSIM file:', fileNESOSIM)
		start = time.time()
		print ('Processing NESOSIM snow depths...')
		dF = cF.gridNESOSIMtoFreeboard(dF, mapProj, fileNESOSIM, dateStr, outSnowVar='snow_depth_N', outDensityVar='snow_density_N')
		print ('Processed NESOSIM snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)

	#------- Get distributed NESOSIM snow depths -----------

	if (nesosimdisttributed==1):
		start = time.time()
		print ('Processing distributed NESOSIM snow depths...')
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_N', outSnowVar='snow_depth_NPdist', consIterations=11, gridSize=100000)
		print ('Processed distributed NESOSIM snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
		
	#-------Thickness conversion-----------
	print ('Processing thickness conversions...')
	if (warrensnow==1):
		#start = time.time()
		#print ('Converting freeboards to thickness...')
		# Convert freeboard to thickness using Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99', rhoi=1)
		#dF.head(3)

	if (modwarrensnow5==1):	
		# Convert freeboard to thickness using modified Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5', rhoi=1)
		#dF.head(3)

	if (nesosim==1):
		#print ('Converting freeboards to thickness 2...')
		# Convert freeboard to thickness using NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_N', snowDensityVar='snow_density_N', outVar='ice_thickness_N', rhoi=1)
		#dF.head(3)

	if (nesosimdisttributed==1):
		#print ('Converting freeboards to thickness 3..')
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdist', rhoi=1)
		#dF.head(3)
		#dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdistUC', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdistUC')
		#dF.head(3)

	print ('Processed thickness conversions')
	#-------Uncertainty calculation-----------
	
	
	# Convert freeboard to thickness using Nathan's distributed NESOSIM data
	#dF = getSnowandsConverttoThickness(dF, snowDepthVar='snowDepthNKdist', snowDensityVar='snowDensityNK', outVar='iceThicknessNKdist')
	#print(dF.head(4))
	#-------Diagnostic plots-----------
	cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+outStr+'W99shot', vars=['freeboard', 'snow_depth_W99', 'snow_density_W99', 'ice_thickness_W99'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod5', 'snowDensityW99mod5', 'iceThicknessW99mod5'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod7', 'snowDensityW99mod7', 'iceThicknessW99mod7'])
	cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+outStr+'Nshot', vars=['freeboard', 'snow_depth_N', 'snow_density_N', 'ice_thickness_N'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NPdistshot', vars=['freeboard', 'snowDepthNPdist', 'snowDensityN', 'iceThicknessNPdist'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'Ndist', vars=['freeboard', 'snowDepthNdist', 'snowDensityN', 'iceThicknessNdist'])
	#plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NKdist', vars=['freeboard', 'snowDepthNKdist', 'snowDensityNK', 'iceThicknessNKdist'])

	#-------Output-----------
	if (savedata==1):
		
		print ('Saving data file to:', dataOutPath+'IS2'+outStr+campaignStr)
		dF.to_pickle(dataOutPath+'IS1'+outStr+campaignStr)
		print ('Saved data')
	if (saveNetCDFX==1):
		
		print ('Saving netCDF (xarray) data file to:', dataOutPath+'IS1'+outStr+campaignStr)
		dFX=xr.Dataset.from_dataframe(dF)
		dFX.to_netcdf(dataOutPath+'IS1'+outStr+campaignStr+'.nc')

if __name__ == '__main__':
	
	iceTypePath='/cooler/scratch1/aapetty/Data/ICETYPE/OSISAF/'
	#snowPath='/cooler/scratch1/aapetty/Data/NESOSIM/OSISAFsig150_ERAI_sf_SICCDR_Rhovariable_IC3_DYN1_WP1_LL1_WPF5.8e-07_WPT5_LLF2.9e-07-100kmnrt3/final/'
	snowPath='/cooler/scratch1/aapetty/Data/NESOSIM/NSIDCv3_ERAI_sf_SICbt_Rhovariable_IC1_DYN1_WP1_LL1_WPF5.8e-07_WPT5_LLF2.9e-07-100kmBSThresh/'
	dataOutPathM='/cooler/scratch1/aapetty/DataOutput/IS1/'
	figPathM='/cooler/scratch1/aapetty/Figures/IS1/'
	ancDataPath='../AncData/'
	
	#global campaignStr
	campaignStr='FM08'
	runStr='run3'

	icesatFreeboardPath='/cooler/scratch1/aapetty/Data/ICESAT/freeboard/'

	#global figPath, dataOutPath
	figPath=figPathM+runStr+'/'+campaignStr+'/'
	dataOutPath=dataOutPathM+runStr+'/'+campaignStr+'/raw/'
	if not os.path.exists(dataOutPath):
		os.makedirs(dataOutPath)
	if not os.path.exists(figPath):
		os.makedirs(figPath)
	
	freeboardFiles = glob(icesatFreeboardPath+'/netcdf/'+campaignStr+'/*.nc')
	print('IS1 campaign season:',campaignStr)
	print('Number of freeboard files: '+str(size(freeboardFiles)))


	#main(freeboardFiles[0])

	with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
		result1=executor.map(main, freeboardFiles)






