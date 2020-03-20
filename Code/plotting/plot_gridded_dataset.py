
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
#mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=33.92, urcrnrlon=102.34, urcrnrlat=31.37)	
def getIS2(savePathT, outStringT, variable):

    dIS2 = xr.open_dataset(savePathT+'IS2ATL10_'+outStringT+'.nc')

    # Get NESOSIM coordinates (constant with time)
    thicknessIS2 = dIS2[variable].values
    latsIS2 = dIS2.latitude.values
    lonsIS2 = dIS2.longitude.values

    xptsT, yptsT = mapProj(lonsIS2, latsIS2)
    return xptsT, yptsT, thicknessIS2 



concDataPath='/cooler/scratch1/aapetty/Data/ICECONC/CDR/monthly/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
iceTypePath='/cooler/scratch1/aapetty/Data/ICETYPE/OSISAF/'

relStr='rel002'
runStr='run12'

figPath='/cooler/scratch1/aapetty/Figures/IS2/'+relStr+'/'+runStr+'/Maps/'
savePath=baseDataPath+'/'+relStr+'/'+runStr+'/products/'
smoothingWindow=200
resolution=25.
beamStr='bnum1'
snowVar='NPdist'
versionStr='vf'
segment=1


variables=['freeboard', 'snow_depth', 'ice_thickness', 'ice_thickness_unc', 'snow_density','ice_type', 'mean_day_of_month', 'num_binned_days']
cbarLabels=['freeboard (m)', 'snow depth (m)', 'ice thickness (m)', 'uncertainity (m)', r'snow density (kg/m$^3$)','ice type', 'mean day of month', 'num valid days in month']
mStr='Apr'
yearStr='2019'
labelStr=runStr+'-'+beamStr+'-'+mStr+yearStr+snowVar+beamStr+'W'+str(smoothingWindow)+'_'+str(resolution)+'km_seg'+str(segment)+versionStr
# DATE INFO
griddedVar=[]

for variableT in variables:
	
	xptsIS2, yptsIS2,griddedVarT = getIS2(savePath, labelStr, variableT)
	griddedVar.append(griddedVarT)

# mask the number of days less than 1
griddedVar[7]=ma.masked_where(griddedVar[7]<1, griddedVar[7])

cbarx=[0.03, 0.276, 0.525, 0.77, 0.03, 0.276, 0.525, 0.77]
cbary=[0.555, 0.555, 0.555, 0.555, 0.08, 0.08, 0.08, 0.08]
minval=[0, 0, 0, 0, 200, -0.5, 0, 0]
maxval=[0.6, 0.4, 4, 1.5, 350,1.5, 30, 20]
numTicks=[4, 5, 5, 4, 4, 2, 4, 5]
numConts=[10, 10, 10, 10, 10, 3, 10, 10]
cmaps=[cm.YlOrRd, cm.BuPu,  cm.cubehelix_r, cm.RdYlBu_r, cm.viridis, cm.Reds, cm.plasma_r, cm.magma_r]

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(8, 5.2))

for i in range(size(variables)):
	ax=axs.flatten()[i]
		
	print(i)
	sca(ax)
	#im1 = mapProj.contourf(xptsIS2 , yptsIS2, ice_thicknessIS2s[i], levels=np.arange(minval, maxval+0.05, 0.5), cmap=cm.viridis , vmin=minval, vmax=maxval,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = mapProj.contourf(xptsIS2 , yptsIS2, griddedVar[i], levels=np.linspace(minval[i], maxval[i], numConts[i]), cmap=cmaps[i], extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	# lower colorbar bounds
	plt.clim(minval[i],maxval[i])
	mapProj.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	mapProj.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)

	mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=5)
	mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	ax.annotate('('+chr(97+i)+') '+variables[i], xy=(0.03, 0.03),backgroundcolor='w', xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=8, zorder=10)

	cax = fig.add_axes([cbarx[i], cbary[i], 0.2, 0.025])

	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
	cbar.set_label(cbarLabels[i], labelpad=1)
	#cbar.ax.locator_params(nbins=4)

	cbar.set_ticks(np.linspace(minval[i], maxval[i], numTicks[i]) )
	if (i==5):
		cbar.set_ticklabels(['FYI', 'MYI'])
	
#ax=axs.flatten()[7]
#fig.delaxes(ax)


subplots_adjust(bottom=0.08, left=0.01, top = 0.99, right=0.99, wspace=0.02, hspace=0.08)
savefig(figPath+'/icethickness'+labelStr+yearStr+mStr+'dataset.png', dpi=300)
