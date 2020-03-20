""" batch_process_icesat2.py
	
	Processing sea ice thickness with ICESat-2 freeboards
	Initial code written by Alek Petty (01/06/2019)
	
	Input:
		ICESat-2 ATL10 freeboards
		NESOSIM snow depths/densities
		Warren snow depth/density
		OSISAF sea ice type

	Output:
		Ice thickness using various snow depth assumptions in netCDF (xarray) format.

	Python dependencies:
		See below for the relevant module imports. Of note:
		xarray/pandas
		netCDF4
		matplotlib
		basemap
		common_functions function library
		Optional (for naive parralel processing):
			itertools
			concurrent.futures
	More information on installation is given in the README file.

	Update history:
		01/06/2020: Version 1.
    
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
from glob import glob
import time
import common_functions as cF

from itertools import repeat, product
import concurrent.futures


def main(fileT, beamNum):
	""" Main ICESat-2 processing module 
	
	Convert the ATL10 segment freeboard data to ice thickness using various input assumptions
	
	Args:
		fileT (str): the ATL10 file path
		beamNum (int): ATLAS beam number from 1 to 6 - 1, 3, 5 (2, 4, 6) are always strong (weak) beams.

	"""

	#====================== Configuration ========================
	# NPdist: NESOSIM snow depth/density, piecewise redistributed, 915 kg/m3 ice density
	# GSFC: W99mod5 snow depth/density, 915 kg/m3 ice density
	# CPOM: regional mean W99mod5, 882, 917 kg/m3 ice density
	# JPL: W99mod5 (bus using ice type fraction), 
	# AWI: W99mod5, 882, 917 kg/m3 ice density

	# I'm sure there is a better way of doing this!..For testing I would turn off a lot of these to speed up computation time.
	regionflags=True
	icetype=True
	warrensnow=True
	modwarrensnow5=True
	modwarrensnow7=True
	nesosim=True
	nesosimdisttributed=True
	nesosimdisttributedrho2=True
	nesosimdisttributedrho3=True
	nesosimdisttributedkwok=True
	modwarrendisttributed=True
	modwarren7disttributed=True
	modwarrensnow5distributedrho3=True
	cpomprocessing=True
	nasaprocessing=True
	awiprocessing=True
	uncertaintycalc=True
	# uncertainity calc needs the distributed nesosim (or select another primary thickness variable)
	saveNetCDFX=True
	
	
	# Map projection
	mapProj = Basemap(projection='npstere',boundinglat=55,lon_0=0, resolution='l' , round=False)

	# Get freeboard files
	start = time.time()
	print ('Grabbing freeboard data...')
	dF, fileStr = cF.getATL10FreeboardShotData(fileT, mapProj, beamNum=beamNum)
	
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
	
	beamStr=dF['beamStr'][0]
	print(beamStr)

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
	
	#------- Get modified (%70) Warren snow depths -----------
	if modwarrensnow7:
		start = time.time()
		print ('Processing Warren mod7 snow depths...')
		dF = cF.getWarrenData(dF, outSnowVar='snow_depth_W99mod7', outDensityVar='None', modFactor=0.7)
		print ('Processed Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
	
	#------- Get regional (CPOM) modified (%50) Warren snow depths -----------
	if cpomprocessing:
		start = time.time()
		print ('Processing CPOM Warren snow depths...')
		dF = cF.getWarrenDataCPOM(dF, outSnowVar='snow_depth_W99mod5r', outDensityVar='snow_density_W99r', modFactor=0.5)
		print ('Processed CPOM Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
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

	if modwarrendisttributed:
		start = time.time()
		print ('Processing distributed mod Warren snow depths...')
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_W99mod5', outSnowVar='snow_depth_W99mod5dist', consIterations=11, gridSize=100000)
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_W99mod5r', outSnowVar='snow_depth_W99mod5rdist', consIterations=11, gridSize=100000)
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_W99mod7', outSnowVar='snow_depth_W99mod7dist', consIterations=11, gridSize=100000)
		print ('Processed distributed mod Warren snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)

	if nesosimdisttributed:
		start = time.time()
		print ('Processing distributed NESOSIM snow depths...')
		dF = cF.distributeSnow(dF, inputSnowDepth='snow_depth_N', outSnowVar='snow_depth_NPdist', consIterations=11, gridSize=100000)
		print ('Processed distributed NESOSIM snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)

	if nesosimdisttributedkwok:
		start = time.time()
		print ('Processing kwok distributed NESOSIM snow depths...')
		dF = cF.distributeSnowKwok(dF, inputSnowDepth='snow_depth_N', outSnowVar='snow_depth_Kdist')
		print ('Processed kwok distributed NESOSIM snow depths in '+str(np.round((time.time()-start), 2))+' seconds')
		#dF.head(3)
		
	#-------Thickness conversion-----------
	print ('Processing thickness conversions...')
	if warrensnow:
		# Convert freeboard to thickness using Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99', rhoi=1)

	if modwarrensnow5:	
		# Convert freeboard to thickness using modified Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5', rhoi=1)
		#dF.head(3)

	if modwarrensnow7:
		# Convert freeboard to thickness using modified Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod7', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod7', rhoi=1)
		#dF.head(3)

	if nesosim:
		# Convert freeboard to thickness using NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_N', snowDensityVar='snow_density_N', outVar='ice_thickness_N', rhoi=1)

	if nesosimdisttributed:
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdist', rhoi=1)

	if nesosimdisttributedrho2:
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdistrho2', rhoi=2)

	if nesosimdisttributedrho3:
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N', outVar='ice_thickness_NPdistrho3', rhoi=3)
		
	if nesosimdisttributedkwok:
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_Kdist', snowDensityVar='snow_density_N', outVar='ice_thickness_Kdist', rhoi=1)
		
	if awiprocessing:
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5', snowDensityVar='snow_density_W99', outVar='ice_thickness_awi', rhoi=2)
		
	if nasaprocessing:
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5', snowDensityVar='snow_density_W99', outVar='ice_thickness_nasa', rhoi=1)
	
	if cpomprocessing:
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5r', snowDensityVar='snow_density_W99r', outVar='ice_thickness_cpom', rhoi=2)

	if modwarrendisttributed:
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5dist', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5dist', rhoi=1)
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5rdist', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5rdist', rhoi=1)
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod7dist', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod7dist', rhoi=1)
				
	if modwarrensnow5distributedrho3:	
		# Convert freeboard to thickness using modified Warren data
		dF = cF.getSnowandConverttoThickness(dF, snowDepthVar='snow_depth_W99mod5dist', snowDensityVar='snow_density_W99', outVar='ice_thickness_W99mod5distrho2', rhoi=2)
		#dF.head(3)

	print ('Processed thickness conversions')
	#-------Uncertainty calculation-----------
	
	
	if uncertaintycalc:
		print ('Processing thickness uncertainity...')
		# Convert freeboard to thickness using distributed NESOSIM data
		dF = cF.getThicknessUncertainty(dF, snowDepthVar='snow_depth_NPdist', snowDensityVar='snow_density_N',iceDensityVar='ice_density_1', outVar='ice_thickness_unc')
		#dF.head(3)
		#print ('Converted freeboards to thickness in '+str(np.round((time.time()-start), 2))+' seconds')
		print ('Processed thickness uncertainity')


	#-------Diagnostic plots-----------
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99shot', vars=['freeboard', 'snowDepthW99', 'snowDensityW99', 'iceThicknessW99'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod5', 'snowDensityW99mod5', 'iceThicknessW99mod5'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'W99mod', vars=['freeboard', 'snowDepthW99mod7', 'snowDensityW99mod7', 'iceThicknessW99mod7'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'Nshot', vars=['freeboard', 'snowDepthN', 'snowDensityN', 'iceThicknessN'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NPdistshot', vars=['freeboard', 'snowDepthNPdist', 'snowDensityN', 'iceThicknessNPdist'])
	#cF.plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'Ndist', vars=['freeboard', 'snowDepthNdist', 'snowDensityN', 'iceThicknessNdist'])
	#plotMap4(dF, mapProj, figPath, dateStr+'_F'+str(fileNum)+'NKdist', vars=['freeboard', 'snowDepthNKdist', 'snowDensityNK', 'iceThicknessNKdist'])
	
	print(dF.head(5))
	print(dF.tail(5))
	
	#-------Output-----------
	if (saveNetCDFX==1):
		outStr=fileT.split("/")[-1][:-3]
		print ('Saving netCDF (xarray) data file to:', dataOutPath+'IS2'+outStr+'_bnum'+str(beamNum)+beamStr)
		dFX=xr.Dataset.from_dataframe(dF)

		dFX.to_netcdf(dataOutPath+'IS2'+outStr+'_bnum'+str(beamNum)+beamStr+'.nc')

if __name__ == '__main__':
	
	iceTypePath='/cooler/scratch1/aapetty/Data/ICETYPE/OSISAF/'
	snowPath='/cooler/scratch1/aapetty/Data/NESOSIM/OSISAFsig150_ERAI_sf_SICCDR_Rhovariable_IC4_DYN1_WP1_LL1_WPF5.8e-07_WPT5_LLF2.9e-07-100kmv56/final/'
	dataOutPathM='/cooler/scratch1/aapetty/DataOutput/IS2/'
	figPathM='/cooler/scratch1/aapetty/Figures/IS2/'
	ancDataPath='../AncData/'
	
	releaseStr='rel002'
	runStr='run13'

	ATL10path='/cooler/I2-ASAS/'+releaseStr+'/ATL10-01/'

	global figPath, dataOutPath
	figPath=figPathM+releaseStr+'/'+runStr+'/'
	dataOutPath=dataOutPathM+releaseStr+'/'+runStr+'/raw/'
	if not os.path.exists(dataOutPath):
		os.makedirs(dataOutPath)
	if not os.path.exists(figPath):
		os.makedirs(figPath)
	
	#The -01 is for Northern Hemisphere datafiles (-02 is for Southern Hemisphere)
	ATL10files0 = glob(ATL10path+'/ATL10-01_2018*.h5')
	ATL10files1 = glob(ATL10path+'/ATL10-01_201901*.h5')
	ATL10files2 = glob(ATL10path+'/ATL10-01_201902*.h5')
	ATL10files3 = glob(ATL10path+'/ATL10-01_201903*.h5')
	ATL10files4 = glob(ATL10path+'/ATL10-01_201904*.h5')
	ATL10files=ATL10files0+ATL10files1+ATL10files2+ATL10files3+ATL10files4
	#print('ATL10 files:', ATL10path+'/ATL10-01_'+dateStr+'*.h5')
	print('ATL10 release:',releaseStr)
	print('Processing run:',runStr)
	print('Number of ATL10 files: '+str(np.size(ATL10files)))
	print(ATL10files)


	beamNums=[1, 2, 3, 4, 5, 6]

	# If you don't want to use concurrent futures and just run over one cpu then use this code
	for ATL10file in ATL10files:
		main(ATL10file, beamNums[0])
	
	#with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:

		# args=((campaign, beam) for beam in beams)
		# print(args)
		# itertools.repeat to add a fixed argument
		# Not very elegant but whatever..d
		#result1=executor.map(main, ATL10files, repeat(beamNums[0]))
		#result2=executor.map(main, ATL10files, repeat(beamNums[2]))
		#result3=executor.map(main, ATL10files, repeat(beamNums[4]))
		#result4=executor.map(main, ATL10files, repeat(beamNums[1]))
		#result5=executor.map(main, ATL10files, repeat(beamNums[3]))
		#result6=executor.map(main, ATL10files, repeat(beamNums[5]))
		
		#result4=executor.map(main, ATL10files, repeat(beamNums[3]))
		#result5=executor.map(main, ATL10files, repeat(beamNums[4]))
		#result6=executor.map(main, ATL10files, repeat(beamNums[5]))
	








