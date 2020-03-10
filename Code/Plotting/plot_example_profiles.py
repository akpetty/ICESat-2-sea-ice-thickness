
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

import time
import sys
sys.path.append('../')
import common_functions as cF




releaseStr='rel002'
runStr='run11'


#figPath='../../Figures/'
figPath='/cooler/scratch1/aapetty/Figures/IS2/'+releaseStr+'/'+runStr+'/'
baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
dataPath=baseDataPath+'/'+releaseStr+'/'+runStr+'/raw/'

example=2
#---Kara Sea FYI
if (example==1):
	beam='bnum1'
	dayStr='18'
	monStr='11'
	fNum=-1
	yearStr='2018'
	sectionNum=201
	lat1=78.5
	lat2=78.6
	lon1=87.34
	lon2=87.44

elif (example==2):

#---- Central Arctic MYI

	beam='bnum1'
	dayStr='16'
	monStr='11'
	fNum=-1
	yearStr='2018'
	#sectionNum=501
	lat1=85.7
	lat2=85.8
	lon1=176.1
	lon2=176.4

cols=['freeboard', 'freeboard_sigma', 'ice_type', 'snow_depth_N', 'snow_density_W99', 'snow_depth_NPdist','snow_depth_Kdist', 'snow_density_N', 'ice_thickness_N', 'ice_thickness_NPdist','ice_thickness_Kdist', 'ice_thickness_unc','ice_thickness_uncrandom', 'ice_thickness_uncsys', 'snow_depth_W99mod5', 'snow_depth_W99mod5dist', 'ice_thickness_W99mod5', 'ice_thickness_W99mod5dist', 'ice_thickness_NPdistrho2']
#cols = vars.copy()
cols.append('delta_time')
cols.append('lat')
cols.append('lon')
cols.append('region_flag')

unit='m'
labelStr=runStr+'-'+yearStr+monStr+dayStr+'_'+beam+'bms'

print(cols)

#IS2data = cF.getATL10Shotdata(dataOutPathM, runStr, campaignStr, cols='all', yearStr=yearStr, monStr=monStr, dayStr=dayStr, fNum=fNum, beamStr=beams[beamNum])
IS2data = cF.getProcessedATL10ShotdataNCDF(dataPath, 
		yearStr=yearStr, monStr=monStr, dayStr=dayStr, fNum=fNum, beamStr=beam, vars=cols)

#DO THIS WHERE LON LAT WITHIN BOUNDS OF THE EXAMPLE I HAD BEFORE!!

#IS2data=IS2data.isel(index=xr.DataArray(np.arange(sectionNum*200, (sectionNum*200)+200), dims=['index']))

IS2data=IS2data.where(((IS2data['lat']>lat1)&(IS2data['lat']<lat2)), drop=True)
IS2data=IS2data.where(((IS2data['lon']>lon1)&(IS2data['lon']<lon2)), drop=True)
#IS2data=IS2data.where(((IS2data['lat']>176.1)&(IS2data['lat']<176.2)), drop=True)

#IS2data=IS2data.isel(index=xr.DataArray(np.arange(sectionNum*200, (sectionNum*200)+200), dims=['index']))

#uncertainity=IS2data[['ice_thickness_NPdist', 'ice_thickness_Kdist', 'ice_thickness_NPdistrho2', 'ice_thickness_W99mod5dist']].to_dataframe().std(axis=1)

ice_type=IS2data['ice_type'].values[0]
if (ice_type>0.5):
	iceTypeStr='MYI'
else:
	iceTypeStr='FYI'

region_label=cF.get_region_mask_sect_labels(int(IS2data['region_flag'].values[0]))
lon1=str(np.round(IS2data['lon'].values[0], 2))
lat1=str(np.round(IS2data['lat'].values[0], 2))

mapProj = Basemap(projection='npstere',boundinglat=55,lon_0=0, resolution='l' , round=False)
xtrk, ytrk=mapProj(IS2data['lon'].values, IS2data['lat'].values)
dist=sqrt((xtrk-xtrk[0])**2+(ytrk-ytrk[0])**2)

density=str(int(floor(IS2data['snow_density_N'].values[0])))
densityWarren=str(int(floor(IS2data['snow_density_W99'].values[0])))

titleStr=region_label+' ('+iceTypeStr+') '+lat1+'N, '+lon1+'E,  Snow density (NESOSIM/W99): '+density+'/'+densityWarren+r' kg m$^{-3}$'

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))

sca(axs.flatten()[0])
ax1=gca()
ax1.plot(dist, IS2data['freeboard'].values, 'x-', color='k', alpha=0.8, label='ATL10 '+beam, markersize=3, linewidth=1)
ax1.errorbar(dist, IS2data['freeboard'].values, yerr=IS2data['freeboard_sigma'].values+0.02, fmt='', linestyle='', marker='.', color='k', lw=0.5, capsize=0.5,)

# add panel indicating ice type
ax1.set_ylabel('freeboard (m)')
#ax1.set_yticks(np.arange(0, 1.5, 0.5))
ax1.set_xticklabels([])

sca(axs.flatten()[1]) 
ax2=gca()
#ax2.yaxis.set_label_position("right")
#ax2.yaxis.tick_right()
ax2.plot(dist, IS2data['snow_depth_N'].values, '.-', color='k', alpha=0.9,label='NSIM', markersize=3, linewidth=1.5)
ax2.plot(dist, IS2data['snow_depth_NPdist'].values, '.-', color='b', alpha=0.9, label=r'NSIM$_{rd-pw}$', markersize=3, linewidth=1)
ax2.plot(dist, IS2data['snow_depth_Kdist'].values, '.-', color='c', alpha=0.9,label=r'NSIM$_{rd-sig}$', markersize=3, linewidth=1)
#ax2.plot(dist, IS2data['snow_depth_NPdistUC'].values, '.-', color='c', alpha=0.9, label=r'NESOSIM$_{dist-uc}$', markersize=3, linewidth=0.7)
#ax2.plot(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['snow_depth_W99'].values, '-', color='r', alpha=0.5,label='W99', linewidth=1.5)
ax2.plot(dist, IS2data['snow_depth_W99mod5'].values, '.-', color='m', alpha=0.9,label='W99m5', markersize=3, linewidth=1.5)
ax2.plot(dist, IS2data['snow_depth_W99mod5dist'].values, '.-', color='r', alpha=0.9,label=r'W99m5$_{rd-pw}$', markersize=3, linewidth=1)
#ax2.plot(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['snow_depth_W99mod7'].values, '-', color='b', alpha=0.5,label='mW99 (0.7)', linewidth=1.5)
#ax2.set_yticks(np.arange(0, 0.4, 0.1))
#ax2.set_ylim([0, 0.35])
ax2.set_ylabel('Snow depth (m)')
ax2.set_xticklabels([])

sca(axs.flatten()[2]) 
ax3=gca()
#ax4.yaxis.set_label_position("right")
#ax4.yaxis.tick_right()
#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
#ax3.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_NPdist'].values-IS2data['ice_thickness_NPdist_unc'].values,IS2data['ice_thickness_NPdist'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
#ax3.plot(dist, IS2data['ice_thickness_N'].values, '-',color='c', alpha=0.9, label='NESOSIM', markersize=3, linewidth=1)

ax3.plot(dist, IS2data['ice_thickness_Kdist'].values, '.-', color='c', alpha=0.9,label=r'NSIM$_{rd-sig}$', markersize=3, linewidth=1)
ax3.plot(dist, IS2data['ice_thickness_NPdistrho2'].values, '.--', color='y', alpha=0.9,label=r'NSIM$_{rd-pw,rho2}$', markersize=3, linewidth=1)
#ax3.plot(dist, IS2data['ice_thickness_NPdistrho1'].values, '.-', color='y', alpha=0.9,label=r'NSIM$_{rd-pw,rho3}$', markersize=3, linewidth=1)

#ax3.plot(dist, IS2data['ice_thickness_NPdistUC'].values, '.-',  color='c', alpha=0.9, label=r'NESOSIM$_{dist-uc}$', markersize=3, linewidth=0.7)
#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99'].values, '-',color='r', alpha=0.5, label='W99', linewidth=1.5)
#ax3.plot(dist, IS2data['ice_thickness_W99mod5'].values, '-',color='m', alpha=0.9, label='W99mod5', markersize=3, linewidth=0.7)
ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, '.-',color='r', alpha=0.9, label=r'W99m5$_{rd-pw}$', markersize=3, linewidth=1)
ax3.plot(dist, IS2data['ice_thickness_NPdist'].values, '.-',  color='b', alpha=0.9, label=r'NSIM$_{rd-pw}$', markersize=3, linewidth=1)
#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)

#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
#ax4.set_yticks(np.arange(0, 5, 1))
ax3.set_ylabel('Ice thickness (m)')
ax3.set_xticklabels([])


sca(axs.flatten()[3]) 
ax4=gca()
#ax4.yaxis.set_label_position("right")
#ax4.yaxis.tick_right()
#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
ax4.plot(dist, IS2data['ice_thickness_unc'].values, '.-', alpha=0.9, label='total', color='k', markersize=3, linewidth=1)
ax4.plot(dist, IS2data['ice_thickness_uncrandom'].values, '.-', alpha=0.9, label='random', color='c', markersize=3, linewidth=1)
ax4.plot(dist, IS2data['ice_thickness_uncsys'].values, '.-', alpha=0.9, label='systematic', color='m', markersize=3, linewidth=1)

ax4.set_ylabel('Thickness uncertainity (m)')
ax4.set_xlabel('Distance (m)')

ax1.legend(loc=1, ncol=1, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
ax2.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
ax3.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
ax4.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)

# Common subplot additions
for ax in axs.flatten():
    sca(ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.grid(axis='both', linestyle='--')

ax1.annotate(titleStr, xy=(0.01, 1.01), xycoords='axes fraction', horizontalalignment='middle', verticalalignment='bottom')

subplots_adjust(left = 0.07, right = 0.98, bottom=0.07, top = 0.96, hspace=0.18)
plt.savefig(figPath+'/ts4'+labelStr+runStr+'_F'+str(fNum)+'eg'+str(example)+'shotDatav3.png', dpi=500)
#plt.savefig(figPathM+'/ts4'+campaignStr+'_F'+str(fileNum)+'shotData.pdf')
#plt.savefig(figPathM+'/ts3'+labelStr+runStr+'_F'+str(fNum)+'shotData.png', dpi=500)
#fig.show()

