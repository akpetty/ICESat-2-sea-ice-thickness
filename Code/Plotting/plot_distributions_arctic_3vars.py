import matplotlib, sys
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
from glob import glob
import netCDF4 as nc4
from scipy.interpolate import griddata
from netCDF4 import Dataset
import scipy.stats as st
import matplotlib.ticker as ticker
import time
import sys
sys.path.append('../')
import common_functions as cF
import seaborn as sns
import dask.array as da

cF.reset_matplotlib()

def get_hists(yearStr, monStr, varT, binVals, binWidth, maxValue):
        
        vars=[varT, 'seg_length', 'region_flag', 'ssh_flag']
        #if (monStr=='12'):
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        #else:
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        print(dataOutPath)
        dFbeams = cF.getProcessedATL10ShotdataNCDF(dataOutPath, yearStr=yearStr, ssh_mask=1, monStr=monStr, dayStr=dayStr, vars=vars, fNum=fNum, beamStr=beam)
        print('Got data')
        dFbeams=dFbeams.where(dFbeams[varT]>0.0, drop=True)
        dFbeams=dFbeams.where(~np.isnan(dFbeams[varT]), drop=True)
        dFbeams=dFbeams.where(dFbeams.seg_length>5, drop=True)

        dFbeams=dFbeams.where(dFbeams.seg_length<200, drop=True)
        
        vals=dFbeams[varT][np.isin(dFbeams.region_flag, regions)]
        segs=dFbeams['seg_length'][np.isin(dFbeams.region_flag, regions)]

        weights=segs/segs.sum().values

        #counts[r, m]=vals.count().values
        meansT=(vals*segs).sum().values/segs.sum().values

        h, bins = da.histogram(vals.data, bins=size(binVals)-1, range=[0, maxValue], weights=weights.data)
        #histVals[m]=h
        histValsT=h.compute()

        return histValsT, meansT

releaseStr='rel002'
runStr='run10'

beam='bnum1'
dayStr='*'
fNum=-1 # -1== all files
snowVar='NPdist'
#var='freeboard'
variables=['freeboard', 'snow_depth_'+snowVar, 'ice_thickness_'+snowVar]
varStrs=['freeboard', 'snow depth', 'ice thickness']
unit='m'
labelStr=runStr+'-'+'-'+dayStr+'_'+beam+'bms'

figPath='/cooler/scratch1/aapetty/Figures/IS2/'+releaseStr+'/'+runStr+'/Dists/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
dataOutPath=baseDataPath+releaseStr+'/'+runStr+'/raw/'

region_label='Inner Arctic\n('+snowVar+')'

regions=[10, 11, 12, 13, 15]

monStrs=['10', '11', '12', '01', '02', '03', '04']
yearStrs=['2018', '2018', '2018', '2019', '2019', '2019', '2019']

numBins=40
histVals=ma.masked_all((3, size(monStrs), numBins))
counts=ma.masked_all((3, size(monStrs)))
means=ma.masked_all((3, size(monStrs)))
binVals=[]
binWidths=[]
maxValues=[0.6, 0.5, 4]

for x in range(3):
    #maxValue=np.percentile((dFbeams[var].values), 95)
    
    binValsT=np.linspace(0, maxValues[x], numBins+1)

    binVals.append(binValsT)
    binWidthT=binValsT[1]-binValsT[0]
    binWidths.append(binWidthT)

    for m in range(size(monStrs)):

        monStr=monStrs[m]
        histVals[x, m], means[x, m] = get_hists(yearStrs[m], monStr, variables[x], binValsT, binWidthT, maxValues[x])




monthlabels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
cmap = cm.magma
c_colors = cmap(np.linspace(0, 0.9, size(monStrs)))


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(4.5, 7))

for x in range(3):
    ax=axs[x]
    sca(ax)
    for m in range(size(monStrs)):
        im = plt.plot(binVals[x][0:-1]+binWidths[x]/2., histVals[x, m],linestyle='-', linewidth=1.6, label=monthlabels[int(monStrs[m])-1], color=c_colors[m])
        
        #meanStr=str(np.round(means[x, m], 2))
        ax.axvline(x=means[x, m], linestyle='--', linewidth=1.6, color=c_colors[m])
    
    if (x==0):
        ax.annotate(region_label, xy=(0.98, 0.85), xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')
     
    if (x==2):
        leg = ax.legend(loc=1, ncol=2,handlelength=1.5, labelspacing=0.2 , numpoints=1, frameon=False)
        ax.set_xlim(0, maxValues[x])


    #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))  
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.025))

    #ax.grid(which='major', axis='both', linestyle='--')

    ax.set_xlabel(varStrs[x]+' (m)')
    #if (x==1):
    ax.set_ylabel('Probability density')   

subplots_adjust(bottom=0.07, left=0.15, right=0.97, top=0.98, hspace=0.25)
plt.savefig(figPath+'/beamTest_'+labelStr+snowVar+'RegionShotSeg_'+str(size(monStrs))+'InnerArctic.png', dpi=500)






