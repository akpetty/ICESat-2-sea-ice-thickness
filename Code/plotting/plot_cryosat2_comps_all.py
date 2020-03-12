
import matplotlib, sys
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
from scipy.interpolate import griddata
from netCDF4 import Dataset
import netCDF4 as nc4
import time
import dask.array as da
import sys
import os
from glob import glob
sys.path.append('../')
import common_functions as cF
cF.reset_matplotlib()

def getCS2gsfc(dataPathCS2, yearStr, mStr):

	f = Dataset(dataPathCS2+'/GSFC/'+yearStr+'/RDEFT4_'+yearStr+mStr+'15.nc', 'r')
	thicknessCS = f.variables['sea_ice_thickness'][:]
	thicknessCS=ma.masked_where(thicknessCS<0, thicknessCS)
	thicknessCS=ma.masked_where(thicknessCS>15, thicknessCS)
	thicknessCS=ma.masked_where(np.isnan(thicknessCS), thicknessCS)
	#thicknessCS=thicknessCS.filled(0)
	#thicknessCS=ma.masked_where(region_mask!=1, thicknessCS)

	latsCS = f.variables['lat'][:]
	lonsCS = f.variables['lon'][:]

	xptsT, yptsT = mapProj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS

def getCS2cpom(dataPathCS2, yearStr, mStr):
	print(dataPathCS2+'/CPOM/thk_'+yearStr+'_'+mStr+'*')
	files=glob(dataPathCS2+'/CPOM/thk_'+yearStr+'_'+mStr+'*')
	
	print(files[0])
	f = Dataset(files[0], 'r')
	latsCS = f.variables['latitude'][::5]
	lonsCS = f.variables['longitude'][::5]
	
	thicknessCS = f.variables['thickness'][::5]

	xptsT, yptsT = mapProj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS

def getCS2jpl(dataPathCS2, yearStr, mStr):
	flat = open(dataPathCS2+'JPL/psn25lats_v3.dat', 'rb')
	flon = open(dataPathCS2+'JPL/psn25lons_v3.dat', 'rb')
	lats = fromfile(file=flat, dtype='<i4')/100000.0
	latsCS = reshape(lats, [448, 304])
	lons = fromfile(file=flon, dtype='<i4')/100000.0
	lonsCS = reshape(lons, [448, 304])
	
	lonsCS=lonsCS[135:-113, 53:-51]
	latsCS=latsCS[135:-113, 53:-51]

	mon1=loadtxt(dataPathCS2+'JPL/month_'+mStr+'.txt')
	thicknessCS=np.flip(mon1, axis=0) 

	xptsT, yptsT = mapProj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS


def getCS2awi(dataPathCS2, yearStr, mStr):
	print(dataPathCS2+'/AWI/*'+yearStr+mStr+'*')
	files=glob(dataPathCS2+'/AWI/*'+yearStr+mStr+'*')
	
	print(files[0])
	f = Dataset(files[0], 'r')
	latsCS = f.variables['lat'][:]
	lonsCS = f.variables['lon'][:]
	
	thicknessCS = f.variables['sea_ice_thickness'][:]

	xptsT, yptsT = mapProj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS

#---- CONFIG --------

#---------- Define file paths --------------------


#mapProj= Basemap(projection='npstere',boundinglat=60,lon_0=0, resolution='l', round=True)
#mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=33.92, urcrnrlon=102.34, urcrnrlat=31.37)	
mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=48., urcrnrlon=100, urcrnrlat=55.37)


def main(CS2PRODUCT, month):
	
	dataPathIS2='/cooler/scratch1/aapetty/DataOutput/IS2/'
	dataPathCS2='/cooler/scratch1/aapetty/Data/CS2/'

	if (month<7):
		year=2019
	else:
		year=2018
	#CS2PRODUCT='JPL' #CPOM, GSFC, JPL

	mStr='%02d' % (month)
	monLabel=cF.monLabels(month-1)
	yearStr=str(year)

	if (CS2PRODUCT=='AWI'):
		snowVar='awi'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2awi(dataPathCS2, yearStr, mStr)
	elif (CS2PRODUCT=='CPOM'):
		snowVar='cpom'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2cpom(dataPathCS2, yearStr, mStr)
	elif (CS2PRODUCT=='JPL'):
		snowVar='nasa'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2jpl(dataPathCS2, yearStr, mStr)
	elif (CS2PRODUCT=='GSFC'):
		snowVar='nasa'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2gsfc(dataPathCS2, yearStr, mStr)
	else:
		print('no data')

	#ice_thicknessCS2=ice_thicknessCS2-0.4

	relStr='rel002'
	runStr='run10'

	figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'

	beamStr='bnum1'

	savePath=dataPathIS2+'/'+relStr+'/'+runStr+'/products/'
	smoothingWindow=400
	resolution=25.

	dayStr='*'

	segment=1
	versionStr='vInt8'

	#snowVar='NPdist'
	labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr


	xptsIS2, yptsIS2,_, _,ice_thicknessIS2 = cF.getIS2gridded(savePath, labelStr, mapProj)


	ice_thicknessCS2[where(ice_thicknessCS2<0.25)]=np.nan
	ice_thicknessIS2[where(ice_thicknessIS2<0.25)]=np.nan

	ice_thicknessCS2G=griddata((xptsCS2.flatten(), yptsCS2.flatten()), ice_thicknessCS2.flatten(), (xptsIS2, yptsIS2), method='nearest')
	ice_thicknessCS2G=ma.masked_where(~np.isfinite(ice_thicknessIS2), ice_thicknessCS2G)
	ice_thicknessCS2G=ma.masked_where(~np.isfinite(ice_thicknessCS2G), ice_thicknessCS2G)

	ice_thicknessIS2=ma.masked_where(~np.isfinite(ice_thicknessCS2G), ice_thicknessIS2)
	ice_thicknessIS2=ma.masked_where(~np.isfinite(ice_thicknessIS2), ice_thicknessIS2)

	region_mask, xptsI, yptsI = cF.get_region_mask_sect('../../AncData/', mapProj, xypts_return=1)
	regions=[10, 11, 12, 13, 15]
	ice_thicknessIS2=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS2)
	ice_thicknessCS2G=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessCS2G)


	ice_thicknessIS2M = ice_thicknessIS2.flatten()[ice_thicknessIS2.flatten().mask == False]
	ice_thicknessCS2M = ice_thicknessCS2G.flatten()[ice_thicknessCS2G.flatten().mask == False]





	trend, sig, r_a, intercept = cF.correlateVars(ice_thicknessCS2M, ice_thicknessIS2M)
	rStr = '%.2f' % r_a
	rmse=sqrt(mean((ice_thicknessIS2M-ice_thicknessCS2M)**2))
	rmsStr='%.2f' % rmse


	merr=mean(ice_thicknessCS2M-ice_thicknessIS2M)
	merrStr='%.2f' % merr	

	std=np.std(ice_thicknessIS2M+merr-ice_thicknessCS2M)
	stdStr='%.2f' % std

	minval=0
	maxval=5


	fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(8, 3.2))
	plt.subplots_adjust(bottom=0.14, left=0.01, top = 0.95, right=0.99, hspace=0.22, wspace=0.2)
	ax1=axs.flatten()[0]
	sca(ax1)
	im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2, levels=np.arange(minval, maxval+0.1, 0.25), cmap=cm.cubehelix_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	#im1=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2,
	#        cmap=cm.viridis, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
	plt.clim(-0.3,5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)


	cax1 = fig.add_axes([0.1, 0.12, 0.3, 0.035])
	cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
	cbar.set_label('Sea ice thickness (m)', labelpad=3)
	cbar.set_ticks(np.arange(minval, maxval+0.1, 1))
	#cbar.set_ticks(np.arange(0, vmaxs[0]+0.1, 0.2))
	ax1.annotate('(a) '+monLabel+' '+yearStr+' IS-2 '+CS2PRODUCT, xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

	ax2=axs.flatten()[1]
	sca(ax2)  

	im2 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessCS2G, levels=np.arange(minval, maxval+0.1, 0.25), cmap=cm.cubehelix_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	#im2=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessCS2G,
	#        cmap=cm.viridis, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
	plt.clim(-0.3,5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

	#cax2 = fig.add_axes([0.28, 0.12, 0.2, 0.035])
	#cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
	#cbar2.set_label('CS2 thickness (m)', labelpad=3)
	ax2.annotate('(b) '+CS2PRODUCT+' CS-2', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

	#cbar2.set_ticks(np.arange(0, vmaxs[1]+0.1, 0.1))

	ax3=axs.flatten()[2]
	sca(ax3) 

	im3 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessCS2G-ice_thicknessIS2, levels=np.arange(-2, 2+0.1, 0.25), cmap=cm.RdBu_r, vmin=-2, vmax=2, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	
	#im3=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessCS2G-ice_thicknessIS2,
	#        cmap=cm.RdBu_r, vmin=-2, vmax=2, zorder=2, rasterized=True)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

	cax3 = fig.add_axes([0.53, 0.12, 0.2, 0.035])
	cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
	cbar3.set_label('difference (m)', labelpad=3)
	cbar3.set_ticks(np.arange(-2, 2.1, 1))
	ax3.annotate('(c) CS-2 minus IS-2', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')


	ax4 = axs.flatten()[3]
	sca(ax4)
	plt.scatter(ice_thicknessCS2G.flatten(), ice_thicknessIS2.flatten(), marker='x', color='0.2', s=4, alpha=.3)
	#nbins, binEdges, _=plt.hist(elevation, bins=30, linewidth=1.5, histtype='step', color='k', density=True, label='elevation')
	#histVals=binEdges+(binEdges[1]-binEdges[0])
	plt.plot(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1), 'k', ls='-', alpha=.5)

	plt.plot(np.arange(0, 10, 0.1), trend*np.arange(0, 10, 0.1)+intercept, 'k', ls='--', alpha=.8)

	ax4.set_xlabel('CS2 thickness (m)', labelpad=1)
	ax4.set_ylabel('IS2 thickness (m)', labelpad=1)
	ax4.annotate('(d)\n r: '+rStr+'\nBias: '+merrStr+' m'+'\nSD: '+stdStr+' m', xy=(0.02, 0.98), 
			xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')

	#ax3.annotate('(c) Elevation distribution' , xy=(0., 1.02), xycoords='axes fraction', color='k', horizontalalignment='left', verticalalignment='bottom')
	ax4.set_xlim(0, 5)
	ax4.set_ylim(0, 5)

	fig.savefig(figPath+'/thicknessComp_'+labelStr+runStr+'CS2IS2corr4'+CS2PRODUCT+'3NP.png', dpi=300)


if __name__ == '__main__':
	#months=[11]
	#for month in months:
	#	main('GSFC', month)
	products=[ 'CPOM', 'JPL', 'AWI', 'GSFC']
	for product in products:
		for month in months:
			main(product, month)
	#months=[11, 12, 1, 2, 3]
	#for month in months:
	#	main('JPL', month)
	#months=[11, 12, 1, 2, 3]
	#for month in months:
	#	main('AWI', month)
	#months=[11]
	#or month in months:
		#main('GSFC', month)

