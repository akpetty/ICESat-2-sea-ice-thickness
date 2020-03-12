
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
#mapProj = Basemap(projection='npstere',boundinglat=56,lon_0=0, resolution='l', round=False  )
mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=48., urcrnrlon=100, urcrnrlat=55.37)

baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS1/'
concDataPath='/cooler/scratch1/aapetty/Data/ICECONC/CDR/monthly/'

runStr='run3'

figPath='/cooler/scratch1/aapetty/Figures/IS1/'+runStr

region_mask, xptsI, yptsI = cF.get_region_mask_sect('../../AncData/', mapProj, xypts_return=1)
#region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsN, yptsN), method='nearest', rescale=True)

regions=[10, 11, 12, 13, 15]

# DATE INFO
ice_thicknessIS1s=[]
dateStrs=[]
iceconcs=[]

#----IS1
campaignStrs=['FM03', 'FM04', 'FM05', 'FM06', 'FM08']
runStr='run3'
versionStr='vInt4'
snowVar='NPdist'
smoothingWindow=0
resolution=25.


for campaignStr in campaignStrs:
	
	savePathIS1=baseDataPath+'/'+runStr+'/'+campaignStr+'/products/'

	labelStr=campaignStr+snowVar+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_'+versionStr

	print(labelStr)
	xptsIS1, yptsIS1,lonsIS1, latsIS1,ice_thicknessIS1 = cF.getIS1gridded(savePathIS1, labelStr, mapProj)

	ice_thicknessIS1=ma.masked_where(np.isnan(ice_thicknessIS1), ice_thicknessIS1)
	ice_thicknessIS1=ma.masked_where(latsIS1>85.5, ice_thicknessIS1)
	ice_thicknessIS1=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS1)

	ice_thicknessIS1s.append(ice_thicknessIS1)
	#xptsc, yptsc, iceconcT=cF.get_cdr_conc(concDataPath, mapProj, yearStr, mStr)
	#iceconcs.append(iceconcT)


# add the mean
ice_thicknessIS1s.insert(0, ma.mean(ice_thicknessIS1s, axis=0))
#dateStrs.append('Winter mean')
campaignStrs.insert(0, 'mean')

minval=0
maxval=5

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(6, 5.8))

i=0
for i in range(size(campaignStrs)):
	ax=axs.flatten()[i]
	#if i == 5:
	#	ax.set_visible(False)
	#else:
	sca(ax)
	#im1 = mapProj.contourf(xptsIS1 , yptsIS1, ice_thicknessIS1s[i], levels=np.arange(minval, maxval+0.05, 0.5), cmap=cm.viridis , vmin=minval, vmax=maxval,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = mapProj.contourf(xptsIS1 , yptsIS1, ice_thicknessIS1s[i], levels=np.arange(minval, maxval+0.1, 0.25), cmap=cm.cubehelix_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	# lower colorbar bounds
	plt.clim(-0.3,5)
	mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	ax.annotate('('+chr(97+i)+') '+campaignStrs[i], xy=(0.98, 0.935),xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom', fontsize=8, zorder=10)
	#im11 = mapProj.contour(xptsc , yptsc, iceconcs[i],levels=0.15, colors='k', linewidths=0.8, zorder=5, alpha=1)
	meanP=str(np.around(ma.mean(ice_thicknessIS1s[i]), decimals=2))
	ax.annotate('Mean: '+meanP+' m', xy=(0.98, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')


#ADD COLORBAR TO MAP

cax = fig.add_axes([0.3, 0.08, 0.4, 0.03])

cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label('Sea ice thickness (m)', labelpad=2)

cbar.set_ticks(np.arange(minval, maxval+1, 1))


subplots_adjust(bottom=0.12, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.03)
savefig(figPath+'/icethicknessIS1'+labelStr+'.png', dpi=300)
