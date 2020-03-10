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



releaseStr='rel002'
runStr='run10'

beam='bnum1'
dayStr='*'
fNum=-1 # -1== all files

#var='freeboard'
var='ice_thickness_NPdist'
varStr='Sea ice thickness'
maxValue=4
binWidth=0.1

vars=[var, 'seg_length', 'region_flag', 'ssh_flag']

unit='m'
labelStr=runStr+'-'+'-'+dayStr+'_'+beam+'bms'


figPath='/cooler/scratch1/aapetty/Figures/IS2/'+releaseStr+'/'+runStr+'/Dists/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
dataOutPath=baseDataPath+releaseStr+'/'+runStr+'/raw/'

#dFbeams=dFbeams.reset_index(drop=True)
regions=[[15], [9], [8], [10], [11], [12], [13]]
region_labels=[cF.get_region_mask_sect_labels(region[0]) for region in regions]

regions.append([10, 11, 12, 13, 15])
region_labels.append('Inner Arctic')

sns.set(style="white")

#maxValue=np.percentile((dFbeams[var].values), 95)

binVals=np.arange(0, maxValue + binWidth, binWidth)


monStrs=['10', '11', '12', '01', '02', '03', '04']
yearStrs=['2018', '2018', '2018', '2019', '2019', '2019', '2019']
monthlabels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
cmap = cm.magma
c_colors = cmap(np.linspace(0, 0.9, size(monStrs)))

histVals=ma.masked_all((size(regions), size(monStrs), size(binVals)-1))
counts=ma.masked_all((size(regions), size(monStrs)))
means=ma.masked_all((size(regions), size(monStrs)))

for m in range(size(monStrs)):

    monStr=monStrs[m]
    dFbeams = cF.getProcessedATL10ShotdataNCDF(dataOutPath, yearStr=yearStrs[m], monStr=monStrs[m], dayStr=dayStr, ssh_mask=1, vars=vars, fNum=fNum, beamStr=beam)
    print('Got data')
    dFbeams=dFbeams.where(dFbeams.seg_length>5, drop=True)
    dFbeams=dFbeams.where(dFbeams.seg_length<200, drop=True)
    dFbeams=dFbeams.where(dFbeams[var]>0.0, drop=True)
    dFbeams=dFbeams.where(~np.isnan(dFbeams[var]), drop=True)

    for r in range(size(regions)):

        #vals=dFbeams.where(dFbeams.region_flag==regions[r], drop=True)[var]
        vals=dFbeams[var][np.isin(dFbeams.region_flag, regions[r])]
        segs=dFbeams['seg_length'][np.isin(dFbeams.region_flag, regions[r])]
        #segs=dFbeams.where(dFbeams.region_flag==regions[r], drop=True)['seg_length']

        weights=segs/segs.sum().values

        #counts[r, m]=vals.count().values
        means[r, m]=(vals*segs).sum().values/segs.sum().values

        h, bins = da.histogram(vals.data, bins=size(binVals)-1, range=[0, maxValue], weights=weights.data)
        histVals[r, m]=h.compute()



#fig = figure(figsize=(10, 10.5))
fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(9, 8))

for x in range(size(regions)):
    ax=axs.flatten()[x]
    sca(ax)
    
    for m in range(size(monStrs)):
        im = plt.plot(binVals[0:-1]+binWidth/2., histVals[x, m],linestyle='-', color=c_colors[m])
        if (x==1):
            im = plt.plot(binVals[0:-1]+binWidth/2., histVals[x, m],linestyle='-', label=monthlabels[int(monStrs[m])-1], color=c_colors[m])
            leg = ax.legend(loc=7, ncol=2,handlelength=1.5, labelspacing=0.2 , numpoints=1, frameon=False)
        
        #meanStr=str(np.round(means[x, m], 2))
        ax.axvline(x=means[x, m], linestyle='--', color=c_colors[m])
        
        #ax.annotate('Mean: '+meanStr+' '+unit+'\n# segs:'+str(counts[x, m]), xy=(0.98, 0.75), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')
   

    ax.annotate('('+chr(97+x)+') '+region_labels[x], xy=(0.98, 0.85), xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')
    ax.set_xlim(0, maxValue)
   

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
     
    #ax.set_ylim(0, 0.15)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    
    ax.grid(which='major', axis='both', linestyle='--')

    if (x>5):
        ax.set_xlabel(varStr+' (m)')
    else:
        ax.set_xticklabels([])
    
    if (x%2 == 0):
        ax.set_ylabel('Probability density') 
    else:
        ax.set_yticklabels([])
    
    
    #if (x==0):
    #    ax.annotate(dateStr, xy=(0.98, 0.4), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')
        #ax.annotate('Bin width: '+str(binWidth)+' m', xy=(0.98, 0.3), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')
        
subplots_adjust(bottom=0.07, left=0.075, right=0.99, top=0.98, hspace = 0.08, wspace=0.11)
plt.savefig(figPath+'/beamTest_'+labelStr+var+'RegionsShotSeg_'+str(size(monStrs))+'monthsv2.png', dpi=500)






