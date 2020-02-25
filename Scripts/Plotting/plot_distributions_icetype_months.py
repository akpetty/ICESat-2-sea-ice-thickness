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


releaseStr='rel001'
runStr='run9'

beam='bnum1'
dayStr='20'
fNum=-1 # -1== all files


var='ice_thickness_NPdist'


relStr='rel002'
runStr='run10'

figPath='../../Figures/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
dataPath=baseDataPath+'/'+relStr+'/'+runStr+'/raw/'
concDataPath='/cooler/scratch1/aapetty/Data/ICECONC/CDR/monthly/'


ice_type_labels=['All', 'MYI', 'FYI']


#maxValue=np.percentile((dFbeams[var].values), 95)
maxValue=20
binWidth=0.1
binVals=np.arange(0, maxValue + binWidth, binWidth)
monStrs=['10', '11', '12', '01', '02', '03', '04']
yearStrs=['2018', '2018', '2018', '2019', '2019', '2019', '2019']
monthlabels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']



histVals=ma.masked_all((size(monStrs), size(ice_type_labels), size(binVals)-1))
#vals=ma.masked_all((size(monStrs), size(ice_type_labels), size(binVals)-1))
means=ma.masked_all((size(monStrs), size(ice_type_labels)))

for m in range(size(monStrs)):
    monStr=monStrs[m]
    dFbeams = cF.getProcessedATL10ShotdataNCDF(dataPath, yearStr=yearStrs[m], monStr=monStrs[m], dayStr=dayStr, fNum=fNum, beamStr=beam)
    dFbeams=dFbeams.where(dFbeams.seg_length>5, drop=True)
    dFbeams=dFbeams.where(dFbeams.seg_length<200, drop=True)
    #dFbeams=dFbeams.where(dFbeams[var]>0.01, drop=True)
    dFbeams=dFbeams.where(dFbeams[var]<20, drop=True)
    dFbeams=dFbeams.where(~np.isnan(dFbeams[var]), drop=True)

    print('Got data')
    for x in range(size(ice_type_labels)):
        if (x==0):
            valsRegion=dFbeams[var]
            segs=dFbeams['seg_length']
        elif (x==1):
            valsRegion=dFbeams.where(dFbeams.ice_type==1, drop=True)[var]
            segs=dFbeams.where(dFbeams.ice_type==1, drop=True)['seg_length']
        else:
            valsRegion=dFbeams.where(dFbeams.ice_type==0, drop=True)[var]
            segs=dFbeams.where(dFbeams.ice_type==0, drop=True)['seg_length']

        if (size(valsRegion)<100):
            continue
        
        weights=segs/segs.sum().values

        #counts[r, m]=vals.count().values
        means[m, x]=(valsRegion*segs).sum().values/segs.sum().values

        h, bins = da.histogram(valsRegion.data, bins=size(binVals)-1, range=[0, maxValue], weights=weights.data)
        histVals[m, x]=h.compute()

savetxt(figPath+'/icetypestats.txt', means)


unit='m'
labelStr=runStr+'-'+dayStr+'_'+beam+'bms'

sns.set(style="white")
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 7))
cmap = cm.magma
c_colors = cmap(np.linspace(0, 0.9, size(monStrs)))


for x in range(size(ice_type_labels)):
    ax=axs.flatten()[x]
    sca(ax)
    ax.annotate('('+chr(97+x)+') '+ice_type_labels[x], xy=(0.98, 0.85), xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')
    ax.set_xlim(0, 4)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.set_ylim(0, 0.15)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))

    #ax.grid(which='minor', axis='both')
    ax.grid(which='major', axis='both', linestyle='--')
     #ax.set_ylim(0, 3)
    if (x>1):
        ax.set_xlabel('Sea ice thickness (m)')
    else:
        ax.set_xticklabels([])
    
    #if (x%2 == 0):
    ax.set_ylabel('Probability density') 
    #else:
        #ax.set_yticklabels([])
       

    for m in range(size(monStrs)):
        
        plt.plot(binVals[0:-1]+(binWidth/2.0), histVals[m, x], color=c_colors[m])
        #plt.plot(binVals[0:-1], vals[m, x], ls='-', color=c_colors[m], linewidth=1)

        if (x==0):
            #plt.plot(binVals[0:-1], vals[m, x], ls='-', label=monthlabels[int(monStrs[m])-1], color=c_colors[m], linewidth=1)

            plt.plot(binVals[0:-1]+(binWidth/2.0), histVals[m, x], linewidth=1.5, label=monthlabels[int(monStrs[m])-1], color=c_colors[m])
            leg = ax.legend(loc=7, ncol=2,handlelength=1.5, labelspacing=0.2 , numpoints=1, frameon=False)
         
        ax.axvline(x=means[m, x], linestyle='--', color=c_colors[m])
        
        
#plt.tight_layout()
subplots_adjust(bottom=0.08, left=0.13, right=0.98, top=0.98, hspace = 0.1, wspace=0.1)
plt.savefig(figPath+'/icetypedistribution_'+labelStr+var+'months.png', dpi=500)
#plt.savefig(figPathM+'/freeboardBeamTest_'+campaignStr+'_F'+str(fileNum)+'shotData.png', dpi=500)




