
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
import rasterio
from rasterio.fill import fillnodata
from calendar import monthrange
from itertools import repeat, product
import concurrent.futures
from datetime import datetime
sys.path.append('../')
import common_functions as cF
cF.reset_matplotlib()

#---------- Define file paths --------------------

dataPathIS2='/cooler/scratch1/aapetty/DataOutput/IS2/'
dataPathCS2='/cooler/scratch1/aapetty/Data/CS2/'
figPathM='/cooler/scratch1/aapetty/Figures/'


def plotData(figPathMv, mPS, xpts1d, ypts1d, var1d, xpts2d, ypts2d, var2d, outStr):
	
	minval=0
	maxval=5
	#for plotting mask data below 0.1 m
	#thicknessG=ma.masked_where(np.isnan(thicknessG), thicknessG)
	minval=np.nanpercentile(var2d, 0)
	maxval=np.nanpercentile(var2d, 95)
	print(minval, maxval)
	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4.5))
	ax1=axs.flatten()[0]
	sca(ax1)
	im1=mPS.hexbin(xpts1d, ypts1d, C=var1d, 
	       cmap=cm.viridis, vmin=minval, vmax=maxval, zorder=2, rasterized=True)


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

	im2=mPS.pcolormesh(xpts2d, ypts2d, var2d, 
	       vmin=minval, vmax=maxval, cmap=cm.viridis, zorder=2, rasterized=True)
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

	dsPS.to_netcdf(savePathT+'IS2ATL10_'+outStringT+'grid.nc')


def outputAlongTrackData(savePathT, lats1d, lons1d, vars1d, monString, outStringT, smoothingWindow):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data will be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  

	Output:

    
    """
     
	f = nc4.Dataset(savePathT+'IS2ATL10_'+outStringT+'alongtrackmean'+'.nc','w', format='NETCDF4') #'w' stands for write
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
	f.author = "Alek Petty, Nathan Kurtz, Ron Kwok, Thorsten Markus"
	f.contact = " alek.a.petty@nasa.gov"
	f.description = "ICESat-2 sea ice thickness"
	
	from datetime import datetime
	f.history = "Created " + datetime.today().strftime("%d/%m/%y")
	#f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def outputBinnedData(savePathT, xpts2d, ypts2d, lats2d, lons2d, vars2d, monString, outStringT):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  

	Output:

    
    """

	f = nc4.Dataset(savePathT+'IS2ATL10_'+outStringT+'.nc','w', format='NETCDF4') #'w' stands for write
	#tempgrp = f.createGroup('DragData')
	print ('dimensions:', lons2d.shape[0], lons2d.shape[1])
	f.createDimension('x', lons2d.shape[0])
	f.createDimension('y', lons2d.shape[1])

	longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
	latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  

	xpts = f.createVariable('xgrid', 'f4', ('x', 'y'))
	ypts = f.createVariable('ygrid', 'f4', ('x', 'y'))  
	

	
	
	longitude[:] = np.around(lons2d, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats2d, decimals=4)
	
	xpts[:] = np.around(xpts2d, decimals=4) #The "[:]" at the end of the variable instance is necessary
	ypts[:] = np.around(ypts2d, decimals=4)

	print(vars2d)
	print(' Variable length:', len(vars2d))


	# Assume all three main variables
	iceThickness = f.createVariable('ice_thickness', 'f4', ('x', 'y'))
	iceThicknessUnc = f.createVariable('ice_thickness_unc', 'f4', ('x', 'y'))
	numDaysData = f.createVariable('num_binned_days', 'i4', ('x', 'y'))
	meanDayOfMonth = f.createVariable('mean_day_of_month', 'f4', ('x', 'y'))
	snowDepth = f.createVariable('snow_depth', 'f4', ('x', 'y'))
	snowDensity = f.createVariable('snow_density', 'f4', ('x', 'y'))
	freeboard = f.createVariable('freeboard', 'f4', ('x', 'y'))
	iceType = f.createVariable('ice_type', 'f4', ('x', 'y'))
	#projection = f.createVariable('projection', 'i4')

	iceThickness[:] = np.around(vars2d[0], decimals=4)
	freeboard[:] = np.around(vars2d[1], decimals=4)
	snowDepth[:] = np.around(vars2d[2], decimals=4)
	snowDensity[:] = np.around(vars2d[3], decimals=4)
	iceType[:] = np.around(vars2d[4], decimals=0)
	iceThicknessUnc[:] = np.around(vars2d[5], decimals=4)
	numDaysData[:] = vars2d[6]
	meanDayOfMonth[:] = np.around(vars2d[7], decimals=0)

	longitude.units = 'degrees East'
	longitude.long_name='longitude'

	xgrid.units = 'meters'
	xgrid.long_name='projection grid in x direction'

	ygrid.units = 'meters'
	ygrid.long_name='projection grid in y direction'

	latitude.units = 'degrees North'
	latitude.long_name='latitude'

	iceThickness.description = monString+' mean ice thickness using piecewise redistributed NESOSIM data.'
	iceThickness.units = 'meters'
	iceThickness.long_name ='sea ice thickness'
	#iceThickness.grid_mapping ='projection'

	iceThicknessUnc.description = monString+' mean ice thickness uncertainty'
	iceThicknessUnc.units = 'meters'
	iceThicknessUnc.long_name ='sea ice thickness uncertainty'
	#iceThicknessUnc.grid_mapping ='projection'

	snowDepth.description = monString+' mean snow depth using piecewise redistributed NESOSIM data.'
	snowDepth.units = 'meters'
	snowDepth.long_name = 'snow depth'
	#snowDepth.grid_mapping ='projection'

	snowDensity.description = monString+' mean snow density from NESOSIM.'
	snowDensity.units = 'kilograms per meters cubed'
	snowDensity.long_name = 'snow density'
	#snowDensity.grid_mapping ='projection'

	freeboard.description = monString+' mean freeboard from ATL10'
	freeboard.units = 'meters'
	freeboard.long_name = 'sea ice freeboard'
	#freeboard.grid_mapping ='projection'

	iceType.description = monString+' mean ice type '
	iceType.units = 'ice type flag (0 = first-year ice, 1 = multi-year ice)'
	iceType.long_name = 'sea ice type classification'
	#iceType.grid_mapping ='projection'

	meanDayOfMonth.description = monString+' mean day of the month represented by a given grid-cell based on the date of the along-track data'
	meanDayOfMonth.units = 'day'
	meanDayOfMonth.long_name = 'day of month'
	#meanDayOfMonth.grid_mapping ='projection'

	numDaysData.description = 'Number of valid daily binnned data in the month'
	numDaysData.units = 'days'
	numDaysData.long_name = 'number of valid daily bins'
	#numDaysData.grid_mapping ='projection'

	#projection.proj4text = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs";
	#projection.spatial_ref='"PROJCS[\"NSIDC Sea Ice Polar Stereographic North\",GEOGCS[\"Unspecified datum based upon the Hughes 1980 ellipsoid\",DATUM[\"Not_specified_based_on_Hughes_1980_ellipsoid\",SPHEROID[\"Hughes 1980\",6378273,298.279411123061,AUTHORITY[\"EPSG\",\"7058\"]],AUTHORITY[\"EPSG\",\"6054\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4054\"]],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],PROJECTION[\"Polar_Stereographic\"],PARAMETER[\"latitude_of_origin\",70],PARAMETER[\"central_meridian\",-45],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],AUTHORITY[\"EPSG\",\"3411\"],AXIS[\"X\",UNKNOWN],AXIS[\"Y\",UNKNOWN]]";'
	#projection.grid_mapping_name = "polar_stereographic"
	#projection.latitude_of_projection_origin = 90.0
	#projection.standard_parallel = 70.0
	#projection.straight_vertical_longitude_from_pole = 135.0
	#projection.longitude_of_projection_origin = -45.0
	#projection.scaling_factor = 1.0
	#projection.false_easting = 0.0
	#projection.false_northing = 0.0
	#projection.semimajor_radius = 6378273.0
	#projection.semiminor_radius = 6356889.449
	#projection.units = "meters"

  #Add global attributes
	f.author = "Alek Petty, Nathan Kurtz, Ron Kwok, Tom Neumann, Thorsten Markus"
	f.contact = " alek.a.petty@nasa.gov"
	f.description = "Gridded ICESat-2  "+monString+" sea ice thickness (and ancillary) data, calculated using ATL10 freeboards and NESOSIM snow depth/density estimates. Data are projected onto the NSIDC Polar Stereographic projection (EPSG:3411)."
	f.reference = "Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann (2020). Winter Arctic sea ice thickness from ICESat-2 freeboards. J. Geophys. Res. Oceans."
	
	
	f.history = "Created " + datetime.today().strftime("%d/%m/%y")
	#f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def getMonthBinned(dataPathIS2v, figPathMv, variable, x, labelStr, snowVar, yearStr, monStr, numDays, fNum, beamStr, resolution, segment, smoothingWindow, mPS, xptsG, yptsG, lonsG, latsG, latMin=60):
	
	IS2dataRawMonth=array([])
	xptsRawMonth=array([])
	yptsRawMonth=array([])
	lonsRawMonth=array([])
	latsRawMonth=array([])
	print(variable, x)

	#variable='freeboard'
	#if (x==1):
		# this is the thickness just over the ice fraction
	#	ssh_mask=1
	#else:
	#	ssh_mask=0
	
	ssh_mask=0

	varBinMonth=np.zeros((numDays, xptsG.shape[0], xptsG.shape[1]))
	varBinMonth[:]=np.nan
	segBinMonth=np.zeros((numDays, xptsG.shape[0], xptsG.shape[1]))
	segBinMonth[:]=np.nan
	#print(varBinSegMonth)
	for day in range(numDays):
		dayStr='%02d' %(day+1)
		print('\n\n\nday:',dayStr)
		print('get IS2 data...')
		
		try:
			IS2dataDay = cF.getProcessedATL10ShotdataNCDF(dataPathIS2v, 
				yearStr=yearStr, monStr=monStr, dayStr=dayStr, fNum=fNum, ssh_mask=ssh_mask, minseg=4, maxseg=200, beamStr=beamStr, vars=['seg_length', 'lat', 'lon', 'ssh_flag', variable], smoothingWindow=smoothingWindow)
			print('Got IS2 data!')
		except:
			continue
		#print('seg range 1:', np.amin(IS2dataDay['seg_length'].values), np.max(IS2dataDay['seg_length'].values))

		# Replace negative values with 0 for regular gridding, but remove for the 'ice effective' gridding.
		#if (x==1):
		#	IS2dataDay=IS2dataDay.where(IS2dataDay[variable]>0, drop=True)
		#else:
		#	IS2dataDay=IS2dataDay.where(IS2dataDay[variable]>0, 0)
		#if (x==1):
			#IS2dataDay=IS2dataDay.where(IS2dataDay[variable]>0.03, drop=True)

		print('Var shape:', IS2dataDay[variable].count)
		IS2dataDay=IS2dataDay.where(~np.isnan(IS2dataDay[variable]), drop=True)

		print('Var shape after nan removal:', IS2dataDay[variable].count)
		
		print('seg range:', np.amin(IS2dataDay['seg_length'].values), np.max(IS2dataDay['seg_length'].values))
		print('val range:', np.amin(IS2dataDay[variable].values), np.amax(IS2dataDay[variable].values))
		

		start = time.time()
		print('Bin data...')
		xptsRaw, yptsRaw=mPS(IS2dataDay['lon'].values, IS2dataDay['lat'].values)
		segBin, binsSeg = cF.bindataN(xptsRaw,yptsRaw, IS2dataDay['seg_length'].values, xptsG, yptsG, resolution*1000.)
		varBin, bins = cF.bindataSegWeighted(xptsRaw, yptsRaw, IS2dataDay[variable].values, IS2dataDay['seg_length'].values, xptsG, yptsG, resolution*1000.)
		
		end = time.time()
		campaignTime = end-start
		campaignTimestr=str(np.round(campaignTime/(60.), 2))
		print ('Proceesing time (minutes): '+campaignTimestr)
		
		if (variable=='ice_thickness_'+snowVar):
			varBin[where(varBin>30)]=30
			varBin[where(varBin<0.0)]=0
		if (variable=='freeboard'):
			varBin[where(varBin>5)]=5
			varBin[where(varBin<0.0)]=0

		minValsinBin=2
		varBin[where(bins<minValsinBin)]=np.nan
		segBin[where(bins<minValsinBin)]=np.nan

		# limit gridded data to latmin (default of 60N)
		#varBin[where(latsG<latMin)]=np.nan

		print('Plot day of binned data...')
		plotData(figPathMv, mPS, xptsRaw, yptsRaw, IS2dataDay[variable].values, xptsG, yptsG, varBin, labelStr+snowVar+dayStr+variable+str(x)+str(smoothingWindow)+str(resolution)+'km_seg'+str(segment))
		#plotData(figPathMv, mPS, xptsRaw, yptsRaw, IS2dataDay['seg_length'].values, xptsG, yptsG, varBinSeg, labelStr+dayStr+'seg_length'+str(smoothingWindow)+str(resolution)+'km_seg'+str(segment))

		varBinMonth[day]=varBin
		segBinMonth[day]=segBin

		IS2dataRawMonth=concatenate([IS2dataRawMonth, IS2dataDay[variable].values])
		xptsRawMonth=concatenate([xptsRawMonth, xptsRaw])
		yptsRawMonth=concatenate([yptsRawMonth, yptsRaw])
		lonsRawMonth=concatenate([lonsRawMonth, IS2dataDay['lon'].values])
		latsRawMonth=concatenate([latsRawMonth, IS2dataDay['lat'].values])

		IS2dataDay.close()

		
	

	#print(IS2dataRawMonth)
	print('month')
	print(varBinMonth)


	varSegWeightedBinMonthMean=np.nansum(varBinMonth*segBinMonth, axis=0)/np.nansum(segBinMonth, axis=0)
	segBinMonthMean=np.nanmean(segBinMonth, axis=0)

	
	# Number of valid data days across the month for each grid-cell
	countDays=np.nansum(~np.isnan(varBinMonth), axis=0)
	
	# mask where less than a certian number of days (should be low based on the orbit cycle)
	mindays=1
	varSegWeightedBinMonthMean[where(countDays<mindays)]=np.nan
	segBinMonthMean[where(countDays<mindays)]=np.nan

	if (variable=='ice_thickness_'+snowVar):
		varSegWeightedBinMonthMean[where(varSegWeightedBinMonthMean>30)]=30
		varSegWeightedBinMonthMean[where(varSegWeightedBinMonthMean<0)]=0
	if (variable=='freeboard'):
		varSegWeightedBinMonthMean[where(varSegWeightedBinMonthMean>5)]=5
		varSegWeightedBinMonthMean[where(varSegWeightedBinMonthMean<0)]=0

	monthMask = ma.zeros((varSegWeightedBinMonthMean.shape))
	monthMask[np.isfinite(varSegWeightedBinMonthMean)]=1

	varSegWeightedBinMonthMean = rasterio.fill.fillnodata(varSegWeightedBinMonthMean, mask=monthMask, max_search_distance=2, smoothing_iterations=0)
	
	#from astropy.convolution import convolve
	#from astropy.convolution import Gaussian2DKernel
	#kernel = Gaussian2DKernel(varBinMonthMean, x_stddev=1)
	#varBinMonthMean = convolve(varBinMonthMean, kernel)

	varSegWeightedBinMonthMean[where(latsG>88)]=np.nan

	plotData(figPathMv, mPS, xptsRawMonth, yptsRawMonth, IS2dataRawMonth, xptsG, yptsG, varSegWeightedBinMonthMean, labelStr+snowVar+variable+str(smoothingWindow)+str(resolution)+'km_seg'+str(segment)+'_month_')

	#xptsRawMonth.close()
	#yptsRawMonth.close()
	return varSegWeightedBinMonthMean, varBinMonth, segBinMonthMean, xptsRawMonth, yptsRawMonth, lonsRawMonth, latsRawMonth, IS2dataRawMonth, countDays

def getPSNgrid(mPS):
	flat = open(dataPathCS2+'JPL/psn25lats_v3.dat', 'rb')
	flon = open(dataPathCS2+'JPL/psn25lons_v3.dat', 'rb')
	lats = fromfile(file=flat, dtype='<i4')/100000.0
	latsG = reshape(lats, [448, 304])
	lons = fromfile(file=flon, dtype='<i4')/100000.0
	lonsG = reshape(lons, [448, 304])
	
	xptsG, yptsG = mPS(lonsG, latsG)

	return lonsG, latsG, xptsG, yptsG

def main(date, resolution=25., relStr='rel002', runStr='run12', snowVar='NPdist', beamStr='bnum1', day=-1, 
	fNum=-1, smoothingWindow=200, segment=1):
	
	month=date[0]
	year=date[1]
	print(month, year)
	#print('Projection...')
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
	#print(lonsG)
	#print(latsG)
	#print(xptsG)
	#print(latsG)
	xidxs=xptsG[0, :]
	yidxs=yptsG[:, 0]
	
	print(month, year, resolution)

	#if (day==-1):
		# all days in month
	#	dayStr='*'
	#else:
	#	dayStr=str(day)
	
	monStr='%02d' %month
	monLabel=cF.monLabels(month-1)
	yearStr=str(year)

	numDays=monthrange(year, month)[1]
	print(numDays)
	#numDays=2
	
	dayVarMonth=np.empty((numDays, xptsG.shape[0], xptsG.shape[1]))
	dayVarMonth[:]=np.nan
	#dayCountMonth=np.empty((numDays, xptsG.shape[0], xptsG.shape[1]))
	#dayCountMonth[:]=np.nan

	 
	
	labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel
	
	dataPathIS2v=dataPathIS2+'/'+relStr+'/'+runStr+'/raw/'

	savePath=dataPathIS2+'/'+relStr+'/'+runStr+'/products/'
	if not os.path.exists(savePath):
	    os.makedirs(savePath)

	figPathMv=figPathM+'/IS2/'+relStr+'/'+runStr+'/Gridding/'+snowVar+'/'+monStr+'/'
	if not os.path.exists(figPathMv):
	    os.makedirs(figPathMv)

	if (snowVar[0:3]=='W99'):
		snowDensityVar='W99'
	else:
		snowDensityVar='N'

	#variables=['ice_thickness_'+snowVar,'ice_thickness_'+snowVar, 'freeboard', 'snow_depth_'+snowVar, 'ice_thickness_uncsys']
	variables=['ice_thickness_'+snowVar, 'freeboard', 'snow_depth_'+snowVar, 'snow_density_'+snowDensityVar, 'ice_type', 'ice_thickness_uncsys']
	#variables=['ice_thickness_'+snowVar, 'snow_depth_'+snowVar, 'snow_density_'+snowDensityVar, 'freeboard', 'ice_type']
	versionStr='v22'

	varsBin=[]
	for x in range(size(variables)):
		variable=variables[x]
		print(variable)
		varBinMonthMean, varBinMonth, segBinMonthMean, xptsRawMonth, yptsRawMonth, lonsRawMonth, latsRawMonth, IS2dataRawMonth, countDays =getMonthBinned(dataPathIS2v, figPathMv, variable, x, labelStr, snowVar, yearStr, monStr, numDays, fNum, beamStr, resolution, segment, smoothingWindow, mPS,  xptsG, yptsG,lonsG, latsG)
		
		varsBin.append(varBinMonthMean)
		outputAlongTrackData(savePath, latsRawMonth, lonsRawMonth, IS2dataRawMonth, monStr, labelStr+variable+beamStr+'W'+str(smoothingWindow)+versionStr, smoothingWindow)
			
		if (x==size(variables)-1):
			
			#varsBin.append(countDays)
			
			# Get mean day
			#dayVarMonth=np.copy(varBinMonth)
			for day in range(numDays):
				#dayCountMonth[day][~np.isnan(varBinMonth[day])]=1.
				#dayCountMonth[day][np.isnan(varBinMonth[day])]=0.
				dayVarMonth[day][~np.isnan(varBinMonth[day])]=day+1

			#varCount=np.nansum(dayCountMonth, axis=0)
			varDay=np.nanmean(dayVarMonth, axis=0)
			plotData(figPathMv, mPS, xptsRawMonth, yptsRawMonth, IS2dataRawMonth, xptsG, yptsG, countDays, labelStr+snowVar+'num_days_data'+str(smoothingWindow)+str(resolution)+'km_seg'+str(segment)+'_month_')
			plotData(figPathMv, mPS, xptsRawMonth, yptsRawMonth, IS2dataRawMonth, xptsG, yptsG, varDay, labelStr+snowVar+'mean_day'+str(smoothingWindow)+str(resolution)+'km_seg'+str(segment)+'_month_')

			#varDay[np.isnan(varBinMonthMean)]=np.nan
			#varCount[np.isnan(varBinMonthMean)]=np.nan
			
			#print('varCount:', varCount)
			print('varDay:', varDay)
			varsBin.append(countDays)
			varsBin.append(varDay)

	#varsBin.append(varDay)

	outputBinnedData(savePath, xptsG, yptsG, latsG, lonsG, varsBin, monLabel, labelStr+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr)


#-- run main program
# Basically allows for main to be called from another script more easily.
if __name__ == '__main__':
	#main(month=10, year=2018)
	#main((11, 2018))
	dates=[(11, 2018), (12, 2018), (1, 2019), (2, 2019), (3, 2019), (4, 2019)]
	#dates=[(11, 2018)]
	with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
		result1=executor.map(main, dates)
	
	#main((12, 2018))
	#main(11, 2018, 25.)
	#main(12, 2018, 25.)
	#main(1, 2019, 25.)
	#main(2, 2019, 25.)
	#main(3, 2019, 25.)
	#main(12, 2018, 25.)
	#main(1, 2019, 25.)
	#main(2, 2019, 25.)

