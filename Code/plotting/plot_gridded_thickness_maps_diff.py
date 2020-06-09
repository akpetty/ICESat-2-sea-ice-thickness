
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

def getIS2(savePathT, outStringT):

    dIS2 = xr.open_dataset(savePathT+'IS2ATL10_'+outStringT+'.nc')

    # Get NESOSIM coordinates (constant with time)
    thicknessIS2 = array(dIS2.ice_thickness)
    latsIS2 = array(dIS2.latitude)
    lonsIS2 = array(dIS2.longitude)

    xptsT, yptsT = mapProj(lonsIS2, latsIS2)
    return xptsT, yptsT, thicknessIS2 


relStr='rel002'
runStr='run12'

baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'
savePath=baseDataPath+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=200
resolution=25.
beamStr='bnum1'
dayStr='*'
#dayStr=str(day)
snowVar='NPdist'
segment=1
versionStr='vf'


years=[2018, 2018, 2019, 2019, 2019, 2019]
months=[11, 12, 1, 2, 3, 4]
# DATE INFO
ice_thicknessIS2s=[]
dateStrs=[]
for x in range(size(years)):
	mStr='%02d' % (months[x])
	monLabel=cF.monLabels(months[x]-1)
	yearStr=str(years[x])
	dateStrs.append(monLabel+' '+yearStr)

	labelStr=runStr+'-'+beamStr+'-'+monLabel+yearStr+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr

	print(labelStr)
	xptsIS2, yptsIS2,ice_thicknessIS2 = getIS2(savePath, labelStr)
	ice_thicknessIS2s.append(ice_thicknessIS2)

minval=-1
maxval=1

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(7.5, 5.8))

i=0
for i in range(size(years)):
	ax=axs.flatten()[i]
	if i == 5:
		ax.set_visible(False)
	else:
		sca(ax)
		im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2s[i+1]-ice_thicknessIS2s[i], levels=np.arange(minval, maxval+0.1, 0.125), cmap=cm.RdBu_r, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
		#im1 = mapProj.pcolormesh(xptsIS2 , yptsIS2, ice_thicknessIS2s[i+1]-ice_thicknessIS2s[i], cmap=cm.RdBu_r , vmin=minval, vmax=maxval, shading='gouraud', zorder=4, rasterized=True)

		mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
		mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

		mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
		mapProj.drawcoastlines(linewidth=0.25, zorder=5)
		ax.annotate('('+chr(97+i)+') '+dateStrs[i]+' to '+dateStrs[i+1], xy=(0.02, 0.92),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=8, zorder=10)

	#m.contour(xpts , ypts, region_maskAO, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskCA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskPS, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskNA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)


#ADD COLORBAR TO MAP

cax = fig.add_axes([0.3, 0.08, 0.4, 0.03])

cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label('thickness change (m)', labelpad=2)

cbar.set_ticks(np.arange(minval, maxval+0.1, 0.5))


subplots_adjust(bottom=0.12, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.03)
savefig(figPath+'/icethickness'+labelStr+str(size(years))+'diff.png', dpi=300)
