
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



beamStr='bnum1'
snowVar='NPdist'


#----IS2


is2DataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
is1DataPath='/cooler/scratch1/aapetty/DataOutput/IS1/'
figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'
savePathIS2=is2DataPath+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=400
resolution=25.

year=2019
month=2
mStr='%02d' % (month)
monLabel=cF.monLabels(month-1)
yearStr=str(year)

dayStr='*'

segment=1
versionStrIS2='vInt8'
campaignStrIS2='FM19'

labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStrIS2



xptsIS2, yptsIS2,_, _,ice_thicknessIS2m1 = cF.getIS2gridded(savePathIS2, labelStr, mapProj, poleHole=85.5)

month=3
monLabel=cF.monLabels(month-1)
labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStrIS2

_, _,lonsIS2, latsIS2,ice_thicknessIS2m2 = cF.getIS2gridded(savePathIS2, labelStr, mapProj, poleHole=85.5)


ice_thicknessIS2=ma.mean([ice_thicknessIS2m1, ice_thicknessIS2m2], axis=0)



#----IS1
campaignStr='FM08'
runStr='run3'
versionStr='vInt4'
#snowVar='N'
savePathIS1=is1DataPath+'/'+runStr+'/'+campaignStr+'/products/'
smoothingWindow=0
labelStr=campaignStr+snowVar+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_'+versionStr

_, _,lonsIS1, latsIS1,ice_thicknessIS1 = cF.getIS1gridded(savePathIS1, labelStr, mapProj)

#IS2mask = ma.zeros((ice_thicknessIS2.shape))
#IS2mask[np.isfinite(ice_thicknessIS2)]=1

#IS1mask = ma.zeros((ice_thicknessIS1.shape))
#IS1mask[np.isfinite(ice_thicknessIS1)]=1


#ice_thicknessIS1 = rasterio.fill.fillnodata(ice_thicknessIS1, mask=IS1mask, max_search_distance=2, smoothing_iterations=0)
#ice_thicknessIS2 = rasterio.fill.fillnodata(ice_thicknessIS2, mask=IS2mask, max_search_distance=2, smoothing_iterations=3)

#ice_thicknessIS2zero=np.copy(ice_thicknessIS2)
ice_thicknessIS1=ma.masked_where(np.isnan(ice_thicknessIS1), ice_thicknessIS1)
ice_thicknessIS2=ma.masked_where(np.isnan(ice_thicknessIS2), ice_thicknessIS2)
#ice_thicknessIS2[where(np.isnan(ice_thicknessIS2))]=0

region_mask, xptsI, yptsI = cF.get_region_mask_sect('../../AncData/', mapProj, xypts_return=1)
#region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsN, yptsN), method='nearest', rescale=True)

regions=[10, 11, 12, 13, 15]
ice_thicknessIS1=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS1)
ice_thicknessIS2=ma.masked_where(~np.isin(region_mask, regions), ice_thicknessIS2)


ice_thicknessIS1=ma.masked_where(latsIS1>85.5, ice_thicknessIS1)
ice_thicknessIS2=ma.masked_where(latsIS2>85.5, ice_thicknessIS2)

vmin=0
vmax=5

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(4, 3.5))
plt.subplots_adjust(bottom=0.14, left=0.01, top = 0.95, right=0.99, hspace=0.22, wspace=0.1)

ax2=axs.flatten()[0]
sca(ax2)  

im2 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2, levels=np.arange(vmin, vmax+0.1, 0.25), cmap=cm.cubehelix_r, vmin=vmin, vmax=vmax, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	

#im2=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2,
#        cmap=cm.cubehelix_r, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
mapProj.drawcoastlines(linewidth=0.25, zorder=5)
mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

#cax2 = fig.add_axes([0.28, 0.12, 0.2, 0.035])
#cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
#cbar2.set_label('CS2 thickness (m)', labelpad=3)
ax2.annotate('(g) IS-2 '+campaignStrIS2, xy=(0.98, 0.935), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

meanP=str(np.around(ma.mean(ice_thicknessIS2), decimals=2))
ax2.annotate('Mean:'+meanP+' m', xy=(0.98, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')


cax1 = fig.add_axes([0.1, 0.12, 0.3, 0.035])
cbar=fig.colorbar(im2, cax=cax1, orientation='horizontal',extend='both')
cbar.set_label('Sea ice thickness (m)', labelpad=3)
cbar.set_ticks(np.arange(vmin, vmax+1, 1))
#cbar2.set_ticks(np.arange(0, vmaxs[1]+0.1, 0.1))

ax3=axs.flatten()[1]
sca(ax3) 

im3 = mapProj.contourf(xptsIS2 , yptsIS2, ma.masked_where(abs(ice_thicknessIS2-ice_thicknessIS1)<0.001, ice_thicknessIS2-ice_thicknessIS1), 
	levels=np.arange(-2, 2+0.1, 0.25), cmap=cm.RdBu_r, vmin=-2, vmax=2, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	

#im3=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2zero-ice_thicknessIS1,
#        cmap=cm.RdBu_r, vmin=-2, vmax=2, zorder=2, rasterized=True)
mapProj.drawcoastlines(linewidth=0.25, zorder=5)
mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

cax3 = fig.add_axes([0.63, 0.12, 0.3, 0.035])
cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
cbar3.set_label('difference (m)', labelpad=3)
cbar3.set_ticks(np.arange(-2, 2.1, 1))
ax3.annotate('(h) '+campaignStrIS2+' - FM08', xy=(0.98, 0.935), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')

meanP=str(np.around(ma.mean(ice_thicknessIS2)-ma.mean(ice_thicknessIS1), decimals=2))
meanP2=str(np.around(ma.mean((ice_thicknessIS2)-ma.mean(ice_thicknessIS1))/ma.mean(ice_thicknessIS1)*100., decimals=0))
ax3.annotate('difference:'+meanP+' m ('+meanP2+'%)', xy=(0.98, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right',color='k')


fig.savefig(figPath+'/thicknessComp_'+runStr+campaignStr+campaignStrIS2+'IS2IS2'+versionStr+versionStrIS2+'v2.png', dpi=300)


