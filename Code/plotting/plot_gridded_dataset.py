
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

def getIS2(savePathT, outStringT, variable):

    dIS2 = xr.open_dataset(savePathT+'IS2ATL10_'+outStringT+'.nc')

    # Get NESOSIM coordinates (constant with time)
    thicknessIS2 = array(dIS2[variable])
    latsIS2 = array(dIS2.latitude)
    lonsIS2 = array(dIS2.longitude)

    xptsT, yptsT = mapProj(lonsIS2, latsIS2)
    return xptsT, yptsT, thicknessIS2 



concDataPath='/cooler/scratch1/aapetty/Data/ICECONC/CDR/monthly/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
iceTypePath='/cooler/scratch1/aapetty/Data/ICETYPE/OSISAF/'

relStr='rel002'
runStr='run11'

figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'
savePath=baseDataPath+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=400
resolution=25.
beamStr='bnum1'
dayStr='*'
#dayStr=str(day)
snowVar='NPdist'
versionStr='vInt9'
segment=1


variables=['freeboard', 'snow_depth', 'ice_thickness', 'ice_thickness_unc', 'snow_density','ice_type', 'num_binned_days', 'mean_day_of_month']
cbarLabels=['freeboard (m)', 'snow depth (m)', 'ice thickness (m)', 'uncertainity (m)', 'snow density (kg/m3)','ice type', 'num days in bin', 'mean day of month']
mStr='Apr'
yearStr='2019'
labelStr=runStr+'-'+beamStr+'-'+yearStr+mStr+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr
# DATE INFO
griddedVar=[]

for variableT in variables:
	
	xptsIS2, yptsIS2,griddedVarT = getIS2(savePath, labelStr, variableT)
	griddedVar.append(griddedVarT)
	
cbarx=[0.03, 0.276, 0.525, 0.77, 0.03, 0.276, 0.525, 0.77]
cbary=[0.555, 0.555, 0.555, 0.555, 0.08, 0.08, 0.08, 0.08]
minval=[0, 0, 0, 0, 200, 0, 0, 0]
maxval=[0.6, 0.4, 4, 1.5, 350,1, 20, 31]
numConts=[10, 10, 10, 10, 10, 3, 10, 10]
cmaps=[cm.YlOrRd, cm.BuPu,  cm.cubehelix_r, cm.RdYlBu_r, cm.viridis, cm.Greens, cm.Reds, cm.magma]

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(8, 5.2))

for i in range(size(variables)):
	ax=axs.flatten()[i]
	print(i)
	#if i == 5:
	#	ax.set_visible(False)
	#else:
	sca(ax)
	#im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2s[i], levels=np.arange(minval, maxval+0.05, 0.5), cmap=cm.viridis , vmin=minval, vmax=maxval,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = mapProj.contourf(xptsIS2 , yptsIS2, griddedVar[i], levels=np.linspace(minval[i], maxval[i], numConts[i]), cmap=cmaps[i], extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	# lower colorbar bounds
	plt.clim(minval[i],maxval[i])
	mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	ax.annotate('('+chr(97+i)+') '+variables[i], xy=(0.01, 0.935),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=8, zorder=10)
	#im11 = mapProj.contour(xptsc , yptsc, iceconcs[i],levels=0.15, colors='k', linewidths=0.8, zorder=5, alpha=1)
	#im11 = mapProj.contour(xptsc , yptsc, iceconcs[i],levels=0.5, colors='r', linewidths=0.5, zorder=5, alpha=1)
	#im21 = mapProj.contour(xptst , yptst, icetypes[i],levels=0.5, colors='k', linewidths=0.4, zorder=5, alpha=0.7)

	#m.contour(xpts , ypts, region_maskAO, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskCA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskPS, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	#m.contour(xpts , ypts, region_maskNA, 1, colors='k', linestyles='-',linewidths=0.4, zorder=4)
	cax = fig.add_axes([cbarx[i], cbary[i], 0.2, 0.025])

	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
	cbar.set_label(cbarLabels[i], labelpad=-1)
	#cbar.ax.locator_params(nbins=4)
	cbar.set_ticks([minval[i], maxval[i]])
	#cbar.set_ticklabels(np.linspace(minval[i], maxval[i]+1, numConts[i]))

#ADD COLORBAR TO MAP




subplots_adjust(bottom=0.08, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.08)
savefig(figPath+'/icethickness'+labelStr+yearStr+mStr+'dataset.png', dpi=300)
