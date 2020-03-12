
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
#import rasterio
#from rasterio.fill import fillnodata
sys.path.append('../')
import common_functions as cF
cF.reset_matplotlib()


#mapProj= Basemap(projection='npstere',boundinglat=60,lon_0=0, resolution='l', round=True)
#mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=33.92, urcrnrlon=102.34, urcrnrlat=31.37)	
mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=48., urcrnrlon=100, urcrnrlat=55.37)

relStr='rel002'
runStr='run10'

figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
beamStr='bnum1'
snowVar='NPdist'
versionStr='vInt8'
versionStr2='vInt8'

#----IS2
savePathIS2=baseDataPath+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=400
resolution=25.

year=2018
month=11
mStr='%02d' % (month)
monLabel=cF.monLabels(month-1)
yearStr=str(year)
dayStr='*'
segment=1


labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr


xptsIS2, yptsIS2,lonsIS2, latsIS2,ice_thicknessIS2v1 = cF.getIS2gridded(savePathIS2, labelStr, mapProj, variable='ice_thickness', poleHole=88)


labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr2

_, _,_, _,ice_thicknessIS2v2 = cF.getIS2gridded(savePathIS2, labelStr, mapProj, variable='ice_thickness_i', poleHole=88)



#ice_thicknessIS2zero=np.copy(ice_thicknessIS2)
ice_thicknessIS2v1=ma.masked_where(np.isnan(ice_thicknessIS2v1), ice_thicknessIS2v1)
ice_thicknessIS2v1=ma.masked_where(np.isnan(ice_thicknessIS2v2), ice_thicknessIS2v1)
ice_thicknessIS2v2=ma.masked_where(np.isnan(ice_thicknessIS2v2), ice_thicknessIS2v2)
ice_thicknessIS2v2=ma.masked_where(np.isnan(ice_thicknessIS2v1), ice_thicknessIS2v2)
#ice_thicknessIS2[where(np.isnan(ice_thicknessIS2))]=0

region_mask, xptsI, yptsI = cF.get_region_mask_sect('../../AncData/', mapProj, xypts_return=1)
#region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsN, yptsN), method='nearest', rescale=True)

regions=[10, 11, 12, 13, 15]
ice_thicknessIS2v1=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS2v1)
ice_thicknessIS2v2=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS2v2)


ice_thicknessIS2v1=ma.masked_where(latsIS2>88., ice_thicknessIS2v1)
ice_thicknessIS2v2=ma.masked_where(latsIS2>88., ice_thicknessIS2v2)

vmin=0
vmax=5

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(7, 3.5))
plt.subplots_adjust(bottom=0.16, left=0.01, top = 0.93, right=0.99, hspace=0.22, wspace=0.1)
ax1=axs.flatten()[0]
sca(ax1)
im1=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2v1,
        cmap=cm.cubehelix_r, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
mapProj.drawcoastlines(linewidth=0.25, zorder=5)
mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)


cax1 = fig.add_axes([0.15, 0.11, 0.3, 0.035])
cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
cbar.set_label('Sea ice thickness (m)', labelpad=3)
#cbar.set_ticks(np.arange(0, vmaxs[0]+0.1, 0.2))
ax1.annotate('(a) ICESat-2 '+versionStr, xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

ax2=axs.flatten()[1]
sca(ax2)  

im2=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2v2,
        cmap=cm.cubehelix_r, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
mapProj.drawcoastlines(linewidth=0.25, zorder=5)
mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

#cax2 = fig.add_axes([0.28, 0.12, 0.2, 0.035])
#cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
#cbar2.set_label('CS2 thickness (m)', labelpad=3)
ax2.annotate('(b) ICESat-2 '+versionStr2, xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

#cbar2.set_ticks(np.arange(0, vmaxs[1]+0.1, 0.1))

ax3=axs.flatten()[2]
sca(ax3) 

#im3 = mapProj.contourf(xptsIS2 , yptsIS2, ma.masked_where(abs(ice_thicknessIS2v2-ice_thicknessIS2v1)<0.001, ice_thicknessIS2v2-ice_thicknessIS2v1), 
#	levels=np.arange(-0.5, 0.5+0.1, 0.1), cmap=cm.RdBu_r, vmin=-0.5, vmax=0.5, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	
im3 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2v2-ice_thicknessIS2v1, 
	levels=np.arange(-0.5, 0.5+0.1, 0.1), cmap=cm.RdBu_r, vmin=-0.5, vmax=0.5, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	
#im3=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2zero-ice_thicknessIS1,
#        cmap=cm.RdBu_r, vmin=-2, vmax=2, zorder=2, rasterized=True)
mapProj.drawcoastlines(linewidth=0.25, zorder=5)
mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

cax3 = fig.add_axes([0.7, 0.11, 0.25, 0.035])
cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
cbar3.set_label('difference (m)', labelpad=3)
cbar3.set_ticks(np.arange(-0.5, 0.55, 0.25))
ax3.annotate('(c) difference', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

meanP=str(np.around(ma.mean(ice_thicknessIS2v2)-ma.mean(ice_thicknessIS2v1), decimals=2))
ax3.annotate('mean:'+meanP+' m', xy=(0.98, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')



fig.savefig(figPath+'/thicknessComp_'+runStr+labelStr+'IS2'+versionStr+versionStr2+'i.png', dpi=300)


