
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

#---------- Define file paths --------------------


#mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=269.26, llcrnrlat=45., urcrnrlon=95.34, urcrnrlat=45.37)
mapProj = Basemap(projection='npstere',boundinglat=56,lon_0=0, resolution='l', round=False  )

def getIS2(savePathT, outStringT, var=''):

    dIS2 = xr.open_dataset(savePathT+'IS2ATL10_'+outStringT+'.nc')

    # Get NESOSIM coordinates (constant with time)
    thicknessIS2 = array(dIS2[var])
    latsIS2 = array(dIS2.latitude)
    lonsIS2 = array(dIS2.longitude)

    xptsT, yptsT = mapProj(lonsIS2, latsIS2)
    return xptsT, yptsT, thicknessIS2 


from config import figPathM
from config import dataPathIS2

relStr='rel001'
runStr='run9'


savePath=dataPathIS2+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=50
resolution=25.
beamStr='bnum1'
dayStr='*'
#dayStr=str(day)
snowVar='NPdist'
versionStr='vInt2'
segment=1


figPath=figPathM+'IS2/'+relStr+'/'+runStr+'/Maps/'


years=[2018, 2018, 2019, 2019, 2019, 2019]
months=[11, 12, 1, 2, 3, 4]
# DATE INFO
ice_thicknessIS2s=[]
ice_thickness_uncIS2s=[]
ice_thickness_unc_percentIS2s=[]
dateStrs=[]
for x in range(size(years)):
	mStr='%02d' % (months[x])
	monLabel=cF.monLabels(months[x]-1)
	yearStr=str(years[x])
	dateStrs.append(monLabel+' '+yearStr)

	labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr

	print(labelStr)
	xptsIS2, yptsIS2,ice_thicknessIS2 = getIS2(savePath, labelStr, var='ice_thickness')
	ice_thicknessIS2s.append(ice_thicknessIS2)

	_, _,ice_thickness_uncIS2 = getIS2(savePath, labelStr, var='thickness_uncertainty')
	ice_thickness_uncIS2s.append(ice_thickness_uncIS2)
	ice_thickness_unc_percentIS2s.append(100*ice_thickness_uncIS2/ice_thicknessIS2)


# add the mean
#ice_thickness_uncIS2s.append(np.nanmean(ice_thickness_uncIS2s, axis=0))
#ice_thickness_unc_percentIS2s.append(np.nanmean(ice_thickness_unc_percentIS2s, axis=0))
#dateStrs.append('Winter mean')


minval=0
maxval=1

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(7.5, 5.8))

i=0
for i in range(size(years)):
	ax=axs.flatten()[i]
	#if i == 5:
	#	ax.set_visible(False)
	#else:
	sca(ax)
	#im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2s[i], levels=np.arange(minval, maxval+0.05, 0.5), cmap=cm.viridis , vmin=minval, vmax=maxval,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thickness_uncIS2s[i], levels=np.arange(minval, maxval+0.1, 0.05), cmap=cm.RdYlBu_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	
	mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	ax.annotate('('+chr(97+i)+') '+dateStrs[i], xy=(0.02, 0.02),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=8, zorder=10)

	#m.contour(xpts , ypts, region_maskAO, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskCA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskPS, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskNA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)


#ADD COLORBAR TO MAP

cax = fig.add_axes([0.3, 0.08, 0.4, 0.03])

cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label('Thickness uncertainity (m)', labelpad=2)

cbar.set_ticks(np.arange(minval, maxval+1, 0.2))


subplots_adjust(bottom=0.12, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.03)
savefig(figPath+'/icethickness'+labelStr+str(size(years))+'unc.png', dpi=300)



minval=0
maxval=100
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(7.5, 5.8))

i=0
for i in range(size(years)):
	ax=axs.flatten()[i]
	#if i == 5:
	#	ax.set_visible(False)
	#else:
	sca(ax)
	#im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2s[i], levels=np.arange(minval, maxval+0.05, 0.5), cmap=cm.viridis , vmin=minval, vmax=maxval,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thickness_unc_percentIS2s[i], levels=np.arange(minval, maxval+0.1, 5), cmap=cm.RdYlBu_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	
	mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	ax.annotate('('+chr(97+i+6)+') '+dateStrs[i], xy=(0.02, 0.02),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=8, zorder=10)

	#m.contour(xpts , ypts, region_maskAO, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskCA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskPS, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskNA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)


#ADD COLORBAR TO MAP

cax = fig.add_axes([0.3, 0.08, 0.4, 0.03])

cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label('Thickness uncertainity (%)', labelpad=2)

cbar.set_ticks(np.arange(minval, maxval+1, 10))


subplots_adjust(bottom=0.12, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.03)
savefig(figPath+'/icethickness'+labelStr+str(size(years))+'uncpercent.png', dpi=300)