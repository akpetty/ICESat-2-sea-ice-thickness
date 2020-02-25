
import matplotlib, sys
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
from scipy.interpolate import griddata
import netCDF4 as nc4
import time
import dask.array as da
import sys
import os
from calendar import monthrange
from itertools import repeat, product
import rasterio
from rasterio.fill import fillnodata
import concurrent.futures
sys.path.append('../')
import common_functions as cF
cF.reset_matplotlib()

#---------- Define file paths --------------------

dataPathIS1='/cooler/scratch1/aapetty/DataOutput/IS1/'
dataPathCS2='/cooler/scratch1/aapetty/Data/CS2/'
figPathM='/cooler/scratch1/aapetty/Figures/'

def plotData(figPathMv, mPS, xpts1d, ypts1d, var1d, xpts2d, ypts2d, var2d, outStr):
	
	minval=0
	maxval=5
	#for plotting mask data below 0.1 m
	#thicknessG=ma.masked_where(np.isnan(thicknessG), thicknessG)

	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4.5))
	ax1=axs.flatten()[0]
	sca(ax1)
	im1=mPS.hexbin(xpts1d, ypts1d, C=var1d, 
	       cmap=cm.viridis, zorder=2, rasterized=True)


	mPS.drawcoastlines(linewidth=0.25, zorder=5)
	mPS.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mPS.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mPS.fillcontinents(color='0.9',lake_color='grey', zorder=7)

	#cax = fig.add_axes([0.4, 0.075, 0.2, 0.03])
	cax = fig.add_axes([0.43, 0.2, 0.03, 0.6])
	cbar = colorbar(im1,cax=cax, orientation='vertical', extend='both',use_gridspec=True)
	#cbar.set_label('Ice thickness (m)', labelpad=3)
	#cbar.set_ticks(np.arange(minval, maxval+1, 1))
	
	ax1.annotate('(a) Raw data ', xy=(0.02, 0.95 ),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', zorder=10)

	ax2=axs.flatten()[1]
	sca(ax2)

	im2=mPS.pcolormesh(xpts2d, ypts2d, ma.masked_where(np.isnan(var2d), var2d), 
	       cmap=cm.viridis, zorder=2, rasterized=True)
	mPS.drawcoastlines(linewidth=0.25, zorder=5)
	mPS.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mPS.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mPS.fillcontinents(color='0.9',lake_color='grey', zorder=7)

	cax2 = fig.add_axes([0.91, 0.2, 0.03, 0.6])
	cbar2 = colorbar(im2,cax=cax2, orientation='vertical', extend='both',use_gridspec=True)
	
	ax2.annotate('(b) Gridded data', xy=(0.02, 0.95 ),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', zorder=10)

	subplots_adjust(bottom=0.01, left=0.01, wspace=0.2, hspace=0.02, top=0.99, right=0.9 )
	plt.savefig(figPathMv+'/gridingTest'+outStr+'.png')
	close(fig)

def plotGrouped(mPS, xptsBin, yptsBin, varBin, outStr):
	
	minval=0
	maxval=5
	#for plotting mask data below 0.1 m
	#thicknessG=ma.masked_where(np.isnan(thicknessG), thicknessG)

	fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
	ax1=gca()

	im1=mPS.hexbin(xptsBin, yptsBin, C=varBin, 
	       cmap=cm.viridis, vmin=minval, vmax=maxval, zorder=2, rasterized=True)


	mPS.drawcoastlines(linewidth=0.25, zorder=5)
	mPS.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mPS.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mPS.fillcontinents(color='0.9',lake_color='grey', zorder=7)

	cax = fig.add_axes([0.7, 0.92, 0.18, 0.03])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
	cbar.set_label('Ice thickness (m)', labelpad=3)
	cbar.set_ticks(np.arange(minval, maxval+1, 1))
	#ax1.annotate('(a) Raw data ', xy=(0.02, 0.95 ),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', zorder=10)

	subdplots_adjust(bottom=0.01, left=0.01, wspace=0.01, hspace=0.02, top=0.99, right=0.99 )

	plt.savefig(figPathMv+'/gridingTest'+outStr+'grouped.png')
	close(fig)


def outNCDF(savePathT, lats2d, lons2d, var2d, outStringT, varString='thickness'):
	# Add projection data and x/y points too?

	dsPS = xr.Dataset({varString: (['x', 'y'],  var2d)},
                coords={'lon': (['x', 'y'], lons2d),
                        'lat': (['x', 'y'], lats2d)})

	dsPS.to_netcdf(savePathT+'IS1_'+outStringT+'grid.nc')


def outputAlongTrackData(savePathT, lats1d, lons1d, vars1d, monString, outStringT, smoothingWindow):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data will be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  

	Output:

    
    """
     
	f = nc4.Dataset(savePathT+'IS1_'+outStringT+'alongtrackmean'+'.nc','w', format='NETCDF4') #'w' stands for write
	#tempgrp = f.createGroup('DragData')
	print('dimensions:', len(lons1d))
	f.createDimension('x', len(lons1d))
	

	longitude = f.createVariable('longitude', 'f4', ('x'))
	latitude = f.createVariable('latitude', 'f4', ('x'))  
	
	#dates = f.createDimension('dates', None)
	#startYr = f.createVariable('startYear', 'i2')
	#date_range = f.createVariable('year', 'str')

	longitude.units = 'degrees East'
	latitude.units = 'degrees North'

	
	longitude[:] = np.around(lons1d, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats1d, decimals=4)
	
	#print('Num vars:', vars2d.shape[0])
	#if (len(vars1d)==1):
	iceThickness = f.createVariable('ice_thickness', 'f4', ('x'))

	iceThickness[:] = np.around(vars1d, decimals=4)
	
	iceThickness.description = monString+' along-track ice thickness mean of '+str(smoothingWindow)+'pts'
	
	#else:
		# Assume all three main variables
		#iceThickness = f.createVariable('ice_thickness', 'f4', ('x'))
		#snowDepth = f.createVariable('snow_depth', 'f4', ('x'))
		#snowDensity = f.createVariable('snow_density', 'f4', ('x'))
		#freeboard = f.createVariable('freeboard', 'f4', ('x'))
		#iceType = f.createVariable('ice_type', 'f4', ('x'))

		#iceThickness[:] = np.around(vars1d[0], decimals=4)
		#snowDepth[:] = np.around(vars1d[1], decimals=4)
		#snowDensity[:] = np.around(vars1d[2], decimals=4)
		#freeboard[:] = np.around(vars1d[3], decimals=4)
		#iceType[:] = np.around(vars1d[4], decimals=4)

		#iceThickness.description = monString+' along-track ice thickness, mean of '+str(smoothingWindow)+'pts'
		#snowDepth.description = monString+' along-track snow depth, mean of '+str(smoothingWindow)+'pts'
		#snowDensity.description = monString+' along-track snow density, mean of '+str(smoothingWindow)+'pts'
		#freeboard.description = monString+' along-track freeboard, mean of '+str(smoothingWindow)+'pts'
		#iceType.description = monString+' along-track ice type, mean of '+str(smoothingWindow)+'pts'


	#Add global attributes
	f.author = "Alek Petty, Nathan Kurtz"
	f.contact = " alek.a.petty@nasa.gov"
	f.description = "ICESat sea ice thickness"
	
	from datetime import datetime
	f.history = "Created " + datetime.today().strftime("%d/%m/%y")
	#f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def outputBinnedData(savePathT, lats2d, lons2d, vars2d, monString, outStringT):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  

	Output:

    
    """

	f = nc4.Dataset(savePathT+'IS1_'+outStringT+'.nc','w', format='NETCDF4') #'w' stands for write
	#tempgrp = f.createGroup('DragData')
	print ('dimensions:', lons2d.shape[0], lons2d.shape[1])
	f.createDimension('x', lons2d.shape[0])
	f.createDimension('y', lons2d.shape[1])

	longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
	latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  
	
	#dates = f.createDimension('dates', None)
	#startYr = f.createVariable('startYear', 'i2')
	#date_range = f.createVariable('year', 'str')

	longitude.units = 'degrees East'
	latitude.units = 'degrees North'

	
	longitude[:] = np.around(lons2d, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats2d, decimals=4)
	
	print(' Variable length:', len(vars2d))
	#print('Num vars:', vars2d.shape[0])
	if (len(vars2d)==1):
		iceThickness = f.createVariable('ice_thickness', 'f4', ('x', 'y'))

		iceThickness[:] = np.around(vars2d[0], decimals=4)
		
		iceThickness.description = monString+' mean ice thickness'
	

	else:
		# Assume all three main variables
		iceThickness = f.createVariable('ice_thickness', 'f4', ('x', 'y'))
		numDaysData = f.createVariable('num_binned_days', 'f4', ('x', 'y'))
		meanDayOfMonth = f.createVariable('mean_day_of_month', 'f4', ('x', 'y'))
		snowDepth = f.createVariable('snow_depth', 'f4', ('x', 'y'))
		snowDensity = f.createVariable('snow_density', 'f4', ('x', 'y'))
		freeboard = f.createVariable('freeboard', 'f4', ('x', 'y'))
		iceType = f.createVariable('ice_type', 'f4', ('x', 'y'))
		

		iceThickness[:] = np.around(vars2d[0], decimals=4)
		numDaysData[:] = np.around(vars2d[1], decimals=4)
		meanDayOfMonth[:] = np.around(vars2d[2], decimals=4)
		snowDepth[:] = np.around(vars2d[3], decimals=4)
		snowDensity[:] = np.around(vars2d[4], decimals=4)
		freeboard[:] = np.around(vars2d[5], decimals=4)
		iceType[:] = np.around(vars2d[6], decimals=4)
		

		iceThickness.description = monString+' mean ice thickness'
		numDaysData.description = 'Number of days of valid binned thickenss data in the month'
		snowDepth.description = monString+' mean snow depth'
		snowDensity.description = monString+' mean snow density'
		freeboard.description = monString+' mean freeboard'
		iceType.description = monString+' mean ice type'
		meanDayOfMonth.description = monString+' mean day of the month'



	#Add global attributes
	f.author = "Alek Petty, Nathan Kurtz"
	f.contact = " alek.a.petty@nasa.gov"
	f.description = "ICESat sea ice thickness"
	
	from datetime import datetime
	f.history = "Created " + datetime.today().strftime("%d/%m/%y")
	#f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def getMonthBinned(dataPathIS1v, figPathMv, variable, campaignStr, snowVar, fNum, resolution, smoothingWindow, mPS, xptsG, yptsG, lonsG, latsG):
	


	try:
		IS1dataMonth = cF.getProcessedIS1(dataPathIS1v, campaignStr, vars=['lat', 'lon', variable], fNum=fNum, smoothingWindow=smoothingWindow)
		print('Got IS data!')
	except:
		print('No data')
		pass

	xptsRaw, yptsRaw=mPS(IS1dataMonth['lon'].values, IS1dataMonth['lat'].values)
	varBin, bins, wherebin = cF.bindataN(xptsRaw,yptsRaw, IS1dataMonth[variable].values, xptsG, yptsG, resolution*1000.)

	minValsinBin=2
	varBin[where(bins<minValsinBin)]=np.nan

	print(varBin)


	if (variable=='ice_thickness_'+snowVar):
		varBin[where(varBin>15)]=15
		varBin[where(varBin<0.0)]=0
	if (variable=='freeboard'):
		varBin[where(varBin>3)]=3
		varBin[where(varBin<0.0)]=0

	monthMask = ma.zeros((varBin.shape))
	monthMask[np.isfinite(varBin)]=1

	varBin = rasterio.fill.fillnodata(varBin, mask=monthMask, max_search_distance=2, smoothing_iterations=0)
	varBin[where(latsG>88)]=np.nan

	plotData(figPathMv, mPS, xptsRaw, yptsRaw, IS1dataMonth[variable].values, xptsG, yptsG, varBin, snowVar+variable+str(smoothingWindow)+str(resolution)+'km_'+campaignStr)

	#xptsRawMonth.close()
	#yptsRawMonth.close()
	return varBin, IS1dataMonth['lon'].values, IS1dataMonth['lat'].values, IS1dataMonth

def getPSNgrid(mPS):
	flat = open(dataPathCS2+'JPL/psn25lats_v3.dat', 'rb')
	flon = open(dataPathCS2+'JPL/psn25lons_v3.dat', 'rb')
	lats = fromfile(file=flat, dtype='<i4')/100000.0
	latsG = reshape(lats, [448, 304])
	lons = fromfile(file=flon, dtype='<i4')/100000.0
	lonsG = reshape(lons, [448, 304])
	
	xptsG, yptsG = mPS(lonsG, latsG)

	return lonsG, latsG, xptsG, yptsG

def main(campaignStr='FM08', resolution=25., runStr='run3', snowVar='NPdist', fNum=-1,
	smoothingWindow=0):
	
	# NSIDC projection
	mPS = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=33.92, urcrnrlon=102.34, urcrnrlat=31.37)	
	# Alek versions
	#mPS = Basemap(epsg=3411,resolution='l', llcrnrlon=259.26, llcrnrlat=62., urcrnrlon=90.34, urcrnrlat=62.37)
	#mPS = Basemap(epsg=3411,resolution='l', llcrnrlon=269.26, llcrnrlat=45., urcrnrlon=95.34, urcrnrlat=45.37)

	# Define our 2D grid we want to bin the data onto
	#lonsG1, latsG1, xptsG1, yptsG1, nx, ny=cF.defGrid(mPS, dxRes=resolution*1000.)
	lonsG, latsG, xptsG, yptsG=getPSNgrid(mPS)
	#print(lonsG1)
	#print(latsG1)
	print(lonsG.shape)

	xidxs=xptsG[0, :]
	yidxs=yptsG[:, 0]
	

	dataPathIS1v=dataPathIS1+'/'+runStr+'/'+campaignStr+'/raw/'

	savePath=dataPathIS1+'/'+runStr+'/'+campaignStr+'/products/'
	print(savePath)
	if not os.path.exists(savePath):
	    os.makedirs(savePath)
	
	figPathMv=figPathM+'/'+runStr+'/Gridded/'+campaignStr+'/'
	if not os.path.exists(figPathMv):
	    os.makedirs(figPathMv)

	if (snowVar[0:3]=='W99'):
		snowDensityVar='W99'
	else:
		snowDensityVar='N'

	variables=['ice_thickness_'+snowVar]
	#variables=['ice_thickness_'+snowVar, 'snow_depth_'+snowVar, 'snow_density_'+snowDensityVar, 'freeboard', 'ice_type']
	versionStr='vInt4'
	varsBin=[]
	for variable in variables:
		
		varBinMonthMean, lonsRawMonth, latsRawMonth, IS1dataRawMonth =getMonthBinned(dataPathIS1v, figPathMv, variable, campaignStr, snowVar, fNum, resolution, smoothingWindow, mPS,  xptsG, yptsG, lonsG, latsG)
		varsBin.append(varBinMonthMean)
			

	outputBinnedData(savePath, latsG, lonsG, varsBin, campaignStr, campaignStr+snowVar+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_'+versionStr)



#-- run main program
if __name__ == '__main__':
	#main(month=10, year=2018)
	main(campaignStr='FM03')
	
