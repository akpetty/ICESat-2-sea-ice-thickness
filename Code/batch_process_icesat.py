""" batch_process_icesat.py
	
	Processing sea ice thickness with ICESat freeboards
	Initial code written by Alek Petty (01/06/2020)
	
	Input:
		ICESat freeboards, NESOSIM snow depths/densities, OSISAF ice type

	Output:
		Ice thickness using various snow depth assumptions.

	Python dependencies:
		See below for the relevant module imports. Of note:
		xarray/pandas
		netCDF4
		matplotlib
		basemap

		More information on installation is given in the README file.

	Update history:
		01/06/2020: Version 1.
    
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



def main(freeboardFile):
	""" Main ICESat-1 processing module 
	
	Convert the ICESat freeboard data to ice thickness using various snow loading input estimates
	
	Args:
		fileT (str): the ATL10 file path
	"""

	outStr=freeboardFile.split("/")[-1][:-3]

	regionflags=True
	icetype=False
	warrensnow=True
	modwarrensnow5=False
	nesosim=True
	nesosimdisttributed=True
	saveNetCDFX=True
	
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
	if regionflags:
		print ('Assign NSIDC region mask...')
		dF=cF.assignRegionMask(dF, mapProj, ancDataPath=ancDataPath)
		#print(dF.head(3)) 
		print ('Complete')

	
	# ----- Get colocated ice type data -------
	if icetype:
		start = time.time()
		print ('Processing ice type data...')
		dF = cF.getIceType(dF, iceTypePath, mapProj, res=4, returnRaw=0)
		print ('Processed ice type in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
	
	#------- Get Warren snow depths -----------
	if warrensnow:
		start = time.time()
		print ('Processing Warren snow depths...')
		dF = cF.getWarrenData(dF, outSnowVar='snow_depth_W99', outDensityVar='snow_density_W99')
		print ('Processed Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
	
	#------- Get modified (50%) Warren snow depths -----------
	if modwarrensnow5:
		start = time.time()
		print ('Processing Warren mod5 snow depths...')
		dF = cF.getWarrenData(dF, outSnowVar='snow_depth_W99mod5', outDensityVar='None', modFactor=0.5)
		print ('Processed Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)

	#------- Get NESOSIM snow depths -----------
	if nesosim:
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

	if nesosimdisttributed:
		start = time.time()
		print ('Processing distributed NESOSIM snow depths...')
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_N', outSnowVar='snow_depth_NPdist', consIterations=11, gridSize=100000)
		print ('Processed distributed NESOSIM snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
		
	#-------Thickness conversion-----------
	print ('Processing thickness conversions...')
	if warrensnow:

		# Convert freeboard to thickness using Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99', rhoi=1)

	if modwarrensnow5:	
		# Convert freeboard to thickness using modified Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5', rhoi=1)

	if nesosim:
		# Convert freeboard to thickness using NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_N', snowDensityVar='snow_density_N', outVar='ice_thickness_N', rhoi=1)

	if nesosimdisttributed:
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdist', rhoi=1)


	print ('Processed thickness conversions')


	#-------Diagnostic plots-----------
	cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+outStr+'W99shot', vars=['freeboard', 'snow_depth_W99', 'snow_density_W99', 'ice_thickness_W99'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod5', 'snowDensityW99mod5', 'iceThicknessW99mod5'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod7', 'snowDensityW99mod7', 'iceThicknessW99mod7'])
	cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+outStr+'Nshot', vars=['freeboard', 'snow_depth_N', 'snow_density_N', 'ice_thickness_N'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NPdistshot', vars=['freeboard', 'snowDepthNPdist', 'snowDensityN', 'iceThicknessNPdist'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'Ndist', vars=['freeboard', 'snowDepthNdist', 'snowDensityN', 'iceThicknessNdist'])
	#plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NKdist', vars=['freeboard', 'snowDepthNKdist', 'snowDensityNK', 'iceThicknessNKdist'])

	#-------Output-----------

	if saveNetCDFX:
		
		print ('Saving netCDF (xarray) data file to:', dataOutPath+'IS1'+outStr+campaignStr)
		dFX=xr.Dataset.from_dataframe(dF)
		dFX.to_netcdf(dataOutPath+'IS1'+outStr+campaignStr+'.nc')

if __name__ == '__main__':
	
	# Input relevant file paths here...
	iceTypePath=''
	snowPath=''
	dataOutPathM=''
	figPathM=''
	ancDataPath='../AncData/'
	
	# ICESat campaign
	campaignStr='FM08'

	# Processing run label
	runStr='run3'

	icesatFreeboardPath=''

	figPath=figPathM+runStr+'/'+campaignStr+'/'
	dataOutPath=dataOutPathM+runStr+'/'+campaignStr+'/raw/'
	if not os.path.exists(dataOutPath):
		os.makedirs(dataOutPath)
	if not os.path.exists(figPath):
		os.makedirs(figPath)
	
	freeboardFiles = glob(icesatFreeboardPath+'/netcdf/'+campaignStr+'/*.nc')
	print('IS1 campaign season:',campaignStr)
	print('Number of freeboard files: '+str(size(freeboardFiles)))

	# If you don't want to use concurrent futures and just run over one cpu then uncomment loop over all files
	#main(freeboardFiles[0])
	
	# Use concurrent futures to run code over multiple CPUs (assigning a given ATL10 granule to a given worker)
	with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
		result1=executor.map(main, freeboardFiles)






