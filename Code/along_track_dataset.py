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

import time
import sys
import common_functions as cF
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
cF.reset_matplotlib()

def plot_along_track(IS2data, IS2datashort, outStr):
	
	dist=(IS2data['along_track_distance'].values-IS2data['along_track_distance'].values[0])*0.001

	distr=(IS2datashort['along_track_distance'].values-IS2datashort['along_track_distance'].values[0])*0.001

	mPS = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=33.92, urcrnrlon=102.34, urcrnrlat=31.37)
	#xpts, ypts=mPS(IS2data['lon'].values, IS2data['lat'].values)
	xpts_short, ypts_short=mPS(IS2datashort['lon'].values, IS2datashort['lat'].values)

	#fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(10, 10))
	fig = plt.figure(figsize=(8, 9))
	ax1 = plt.subplot2grid((6,4), (0,0), colspan=3)


	ax1.plot(dist, IS2data['freeboard'].values, 'x', color='k', alpha=0.2, label='ATL10 '+beam, markersize=2, linewidth=1)
	ax1.plot(distr, IS2datashort['freeboard'].values, 'x', color='r', alpha=0.8, label='ATL10 '+beam+' (200 segment rolling mean)', markersize=2, linewidth=1)

	# add panel indicating ice type
	ax1.set_ylabel('freeboard (m)')
	#ax1.set_yticks(np.arange(0, 1.5, 0.5))
	ax1.set_xticklabels([])

	ax2 = plt.subplot2grid((6,4), (1,0), colspan=3)
	#ax2.yaxis.set_label_position("right")
	#ax2.yaxis.tick_right()
	#ax2.plot(dist, IS2data['snow_depth_N'].values, 'x', color='k', alpha=0.2,label='NSIM', markersize=2, linewidth=1.5)
	ax2.plot(dist, IS2data['snow_depth_NPdist'].values, 'x', color='k', alpha=0.2, label=r'NS$_{rd-pw}$', markersize=2, linewidth=1)
	ax2.plot(distr, IS2datashort['snow_depth_NPdist'].values, 'x', color='r', alpha=0.8, label=r'NS$_{rd-pw}$ (200 segment mean)', markersize=2, linewidth=1)

	#ax2.plot(dist, IS2data['snow_depth_W99mod5'].values, 'x', color='m', alpha=0.2,label='W99m5', markersize=2, linewidth=1.5)
	#ax2.plot(dist, IS2data['snow_depth_W99mod5dist'].values, 'x', color='r', alpha=0.2,label=r'W99m5$_{rd-pw}$', markersize=2, linewidth=1)
	#ax2.plot(dist, IS2data['snow_depth_W99mod7dist'].values, 'x', color='y', alpha=0.2,label=r'W99m7$_{rd-pw}$', markersize=2, linewidth=1)
	#ax2.plot(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['snow_depth_W99mod7'].values, '-', color='b', alpha=0.5,label='mW99 (0.7)', linewidth=1.5)
	#ax2.set_yticks(np.arange(0, 0.4, 0.1))
	#ax2.set_ylim([0, 0.35])
	ax2.set_ylabel('Snow depth (m)')
	ax2.set_xticklabels([])

	ax3 = plt.subplot2grid((6,4), (2,0), colspan=3)
	#ax4.yaxis.set_label_position("right")
	#ax4.yaxis.tick_right()
	#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
	#ax3.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_NPdist'].values-IS2data['ice_thickness_NPdist_unc'].values,IS2data['ice_thickness_NPdist'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax3.plot(dist, IS2data['ice_thickness_N'].values, '-',color='c', alpha=0.9, label='NESOSIM', markersize=3, linewidth=1)

	#ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, 'x',color='r', alpha=0.2, label=r'W99m5$_{rd-pw}$', markersize=2, linewidth=1)
	ax3.plot(dist, IS2data['snow_density_N'].values, 'x',  color='k', alpha=0.2, label=r'NS', markersize=2, linewidth=1)
	ax3.plot(distr, IS2datashort['snow_density_N'].values, 'x',  color='r', alpha=0.9, label=r'NS (200 segment mean)', markersize=2, linewidth=1)
	#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)
	#ax3.errorbar(dist, IS2data['ice_thickness_NPdist'].values, yerr=IS2data['ice_thickness_unc'].values, fmt='', linestyle='', marker='.', color='k', lw=0.5, capsize=0.5,)

	#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax4.set_yticks(np.arange(0, 5, 1))
	ax3.set_ylabel(r'Snow density (kg/m$^3$)')
	ax3.set_xticklabels([])

	ax4 = plt.subplot2grid((6,4), (3,0), colspan=3)
	#ax4.yaxis.set_label_position("right")
	#ax4.yaxis.tick_right()
	#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
	#ax3.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_NPdist'].values-IS2data['ice_thickness_NPdist_unc'].values,IS2data['ice_thickness_NPdist'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax3.plot(dist, IS2data['ice_thickness_N'].values, '-',color='c', alpha=0.9, label='NESOSIM', markersize=3, linewidth=1)

	#ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, 'x',color='r', alpha=0.2, label=r'W99m5$_{rd-pw}$', markersize=2, linewidth=1)
	ax4.plot(dist, IS2data['ice_thickness_NPdist'].values, 'x',  color='k', alpha=0.2, label=r'NS$_{rd-pw}$', markersize=2, linewidth=1)
	ax4.plot(distr, IS2datashort['ice_thickness_NPdist'].values, 'x',  color='r', alpha=0.9, label=r'NS$_{rd-pw}$ (200 segment mean)', markersize=2, linewidth=1)
	#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)
	#ax3.errorbar(dist, IS2data['ice_thickness_NPdist'].values, yerr=IS2data['ice_thickness_unc'].values, fmt='', linestyle='', marker='.', color='k', lw=0.5, capsize=0.5,)

	#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax4.set_yticks(np.arange(0, 5, 1))
	ax4.set_ylabel('Ice thickness (m)')
	ax4.set_xticklabels([])

	ax5 = plt.subplot2grid((6,4), (4,0), colspan=3)
	#ax4.yaxis.set_label_position("right")
	#ax4.yaxis.tick_right()
	#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
	#ax3.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_NPdist'].values-IS2data['ice_thickness_NPdist_unc'].values,IS2data['ice_thickness_NPdist'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax3.plot(dist, IS2data['ice_thickness_N'].values, '-',color='c', alpha=0.9, label='NESOSIM', markersize=3, linewidth=1)

	#ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, 'x',color='r', alpha=0.2, label=r'W99m5$_{rd-pw}$', markersize=2, linewidth=1)
	ax5.plot(dist, IS2data['ice_thickness_unc'].values, 'x',  color='k', alpha=0.2, label=r'total uncertainty$', markersize=2, linewidth=1)
	ax5.plot(distr, IS2datashort['ice_thickness_unc'].values, 'x',  color='r', alpha=0.9, label=r'total uncertainty (200 segment mean)', markersize=2, linewidth=1)
	#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)
	#ax3.errorbar(dist, IS2data['ice_thickness_NPdist'].values, yerr=IS2data['ice_thickness_unc'].values, fmt='', linestyle='', marker='.', color='k', lw=0.5, capsize=0.5,)

	#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
	#ax4.set_yticks(np.arange(0, 5, 1))
	ax5.set_ylabel('Uncertainty (m)')
	ax5.set_xticklabels([])


	ax6 = plt.subplot2grid((6,4), (5,0), colspan=3)
	#ax4.yaxis.set_label_position("right")
	#ax4.yaxis.tick_right()
	#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
	ax6.plot(dist, IS2data['seg_length'].values, 'x', alpha=0.4, label='segment length', color='k', markersize=2, linewidth=1)
	ax6.plot(distr, IS2datashort['seg_length'].values, 'x', alpha=0.9, label='segment length (200 segment mean)', color='r', markersize=2, linewidth=1)

	ax6.set_ylabel('Segment length (m)')
	ax6.set_xlabel('Along track distance from x (km)')

	ax1.legend(loc=1, ncol=4, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
	ax2.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
	ax3.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
	ax4.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
	ax5.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
	ax6.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)

	# Common subplot additions
	for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
	    sca(ax)
	    ax.spines['right'].set_visible(False)
	    ax.spines['top'].set_visible(False)
	    ax.xaxis.set_ticks_position('bottom')
	    ax.yaxis.set_ticks_position('left')
	    ax.grid(axis='both', linestyle='--')

	ax1.annotate(outStr, xy=(0.02, 1.2), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

	ax6 = plt.subplot2grid((6,4), (0,3), rowspan=2)
	im0=mPS.plot(xpts_short[0], ypts_short[0], '*', color='k', markersize=7, zorder=5)
	im1=mPS.scatter(xpts_short, ypts_short, c=IS2datashort['ice_thickness_NPdist'].values, s=5,
	        cmap=cm.viridis, vmin=0, vmax=5, zorder=2, rasterized=True)
	mPS.drawcoastlines(linewidth=0.25, zorder=5)
	mPS.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mPS.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	mPS.fillcontinents(color='0.9',lake_color='grey', zorder=3)
	#cax1 = fig.add_axes([0.79, 0.57, 0.18, 0.02])
	#cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
	#cbar.set_label('Sea ice thickness (m)', labelpad=3)

	subplots_adjust(left = 0.07, right = 0.98, bottom=0.06, top = 0.94, hspace=0.18)
	plt.savefig(outPath+outStr+'_along_track.png', dpi=300)



def main(fileT):

	#---- Get Data-----
	print(fileT)
	outStr='IS2seaicethickness_'+((fileT.split('IS2ATL10-'))[1].split('.nc')[0])

	IS2data= xr.open_mfdataset(fileT, concat_dim='index', parallel=True)[cols]
	print(IS2data)
	smoothingWindow=200
	for col in cols:

		IS2dataR=IS2data[col].rolling(index=smoothingWindow, center=True).mean()
		IS2dataRM=IS2dataR[int(smoothingWindow/2):-int(smoothingWindow/2):smoothingWindow]
		if col=='delta_time':
			IS2data_short = IS2dataRM.to_dataset(name = col)
		else:
			IS2data_short[col] = IS2dataRM
	
	print(IS2data_short)

	

	#---- Save Data-----
	
	IS2data['delta_time'].attrs={'units':'seconds', 'long_name':'Delta time'}
	IS2data['along_track_distance'].attrs={'units':'meters', 'long_name':'Along track distance'}
	IS2data['lat'].attrs={'units':'degrees East', 'long_name':'Latitude'}
	IS2data['lon'].attrs={'units':'degrees North', 'long_name':'Longitude'}
	IS2data['ssh_flag'].attrs={'units':'surface flag', 'long_name':'Sea surface flag'}
	IS2data['seg_length'].attrs={'units':'meters', 'long_name':'Segment length'}
	IS2data['freeboard'].attrs={'units':'meters', 'long_name':'Sea ice freeboard'}
	IS2data['snow_depth_N'].attrs={'units':'meters', 'long_name':'NESOSIM snow depth'}
	IS2data['snow_depth_NPdist'].attrs={'units':'meters', 'long_name':'NESOSIM snow depth piecewise redistributed'}
	IS2data['snow_depth_W99mod5'].attrs={'units':'meters', 'long_name':'Modified 50 percent FYI Warren snow depth'}
	IS2data['snow_depth_W99mod7'].attrs={'units':'meters', 'long_name':'Modified 70 percent FYI Warren snow depth'}
	IS2data['snow_depth_W99mod5r'].attrs={'units':'meters', 'long_name':'Modified Warren snow depth cpom'}
	IS2data['snow_depth_W99mod5dist'].attrs={'units':'meters', 'long_name':'Modified 50 percent FYI Warren snow depth piecewise redistributed'}
	IS2data['snow_depth_W99mod7dist'].attrs={'units':'meters', 'long_name':'Modified 70 percent FYI Warren snow depth piecewise redistributed'}
	IS2data['snow_depth_W99mod5rdist'].attrs={'units':'meters', 'long_name':'Modified Warren snow depth cpom piecewise redistributed'}
	IS2data['snow_density_N'].attrs={'units':'kilograms per meter cubed', 'long_name':'NESOSIM snow density'}
	IS2data['snow_density_W99'].attrs={'units':'kilograms per meter cubed', 'long_name':'Warren snow density'}
	IS2data['snow_density_W99r'].attrs={'units':'kilograms per meter cubed', 'long_name':'Warren snow density cpom'}
	IS2data['ice_thickness_NPdist'].attrs={'units':'meters', 'long_name':'Sea ice thickness from redistributed NESOSIM snow loading'}
	IS2data['ice_thickness_W99mod5dist'].attrs={'units':'meters', 'long_name':'Sea ice thickness from redistributed modified 50 percent FYU Warren snow loading'}
	IS2data['ice_thickness_unc'].attrs={'units':'meters', 'long_name':'Total sea ice thickness uncertainty'}
	IS2data['ice_thickness_uncsys'].attrs={'units':'meters', 'long_name':'Systematic sea ice thickness uncertainty'}
	IS2data['ice_thickness_uncrandom'].attrs={'units':'meters', 'long_name':'Random sea ice thickness uncertainty'}
	IS2data['ice_type'].attrs={'units':'type flag', 'long_name':'OSI SAF sea ice type classification'}


	IS2data.attrs={'Title':'Along track sea ice thickness estimates from ICESat-2 and ancillary data including ATL10 freeboards and NESOSIM snow loading.',
	'Author':'Alek Petty, Nathan Kurtz, Ron Kwok, Tom Neumann, Thorsten Markus',
	'Contact':'Alek Petty (alek.a.petty@nasa.gov)',
	'Reference':'Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann (2020). Winter Arctic sea ice thickness from ICESat-2 freeboards. J. Geophys. Res. Oceans.'}
	IS2data.to_netcdf(outPath+outStr+'_along_track.nc')

	IS2data_short['delta_time'].attrs={'units':'seconds', 'long_name':'Delta time'}
	IS2data_short['along_track_distance'].attrs={'units':'meters', 'long_name':'Along track distance'}
	IS2data_short['lat'].attrs={'units':'degrees', 'long_name':'Latitude'}
	IS2data_short['lon'].attrs={'units':'degrees', 'long_name':'Longitude'}
	IS2data_short['ssh_flag'].attrs={'units':'flag', 'long_name':'Sea surface flag'}
	IS2data_short['seg_length'].attrs={'units':'meters', 'long_name':'Segment length'}
	IS2data_short['freeboard'].attrs={'units':'meters', 'long_name':'Sea ice freeboard'}
	IS2data_short['snow_depth_N'].attrs={'units':'meters', 'long_name':'NESOSIM snow depth'}
	IS2data_short['snow_depth_NPdist'].attrs={'units':'meters', 'long_name':'NESOSIM snow depth piecewise redistributed'}
	IS2data_short['snow_depth_W99mod5'].attrs={'units':'meters', 'long_name':'Modified 50 percent FYI Warren snow depth'}
	IS2data_short['snow_depth_W99mod7'].attrs={'units':'meters', 'long_name':'Modified 70 percent FYI Warren snow depth'}
	IS2data_short['snow_depth_W99mod5r'].attrs={'units':'meters', 'long_name':'Modified Warren snow depth cpom'}
	IS2data_short['snow_depth_W99mod5dist'].attrs={'units':'meters', 'long_name':'Modified 50 percent FYI Warren snow depth piecewise redistributed'}
	IS2data_short['snow_depth_W99mod7dist'].attrs={'units':'meters', 'long_name':'Modified 70 percent FYI Warren snow depth piecewise redistributed'}
	IS2data_short['snow_depth_W99mod5rdist'].attrs={'units':'meters', 'long_name':'Modified Warren snow depth cpom piecewise redistributed'}
	IS2data_short['snow_density_N'].attrs={'units':'kilograms per meter cubed', 'long_name':'NESOSIM snow density'}
	IS2data_short['snow_density_W99'].attrs={'units':'kilograms per meter cubed', 'long_name':'Warren snow density'}
	IS2data_short['snow_density_W99r'].attrs={'units':'kilograms per meter cubed', 'long_name':'Warren snow density cpom'}
	IS2data_short['ice_thickness_NPdist'].attrs={'units':'meters', 'long_name':'Sea ice thickness from redistributed NESOSIM snow loading'}
	IS2data_short['ice_thickness_W99mod5dist'].attrs={'units':'meters', 'long_name':'Sea ice thickness from redistributed modified 50 percent FYU Warren snow loading'}
	IS2data_short['ice_thickness_unc'].attrs={'units':'meters', 'long_name':'Total sea ice thickness uncertainty'}
	IS2data_short['ice_thickness_uncsys'].attrs={'units':'meters', 'long_name':'Systematic sea ice thickness uncertainty'}
	IS2data_short['ice_thickness_uncrandom'].attrs={'units':'meters', 'long_name':'Random sea ice thickness uncertainty'}
	IS2data_short['ice_type'].attrs={'units':'type flag', 'long_name':'OSI SAF sea ice type classification'}

	#IS2data_short.attrs={'Title':'Smoothed (200 segment) along track sea ice thickness estimates from ICESat-2 and ancillary data', 'Data summary':'Smoothed (200 segment) along track data generated during the processing of ATL10 freebaords into sea ice thickness. Data generated by Alek Petty (alek.a.petty@nasa.gov).'}
	IS2data_short.attrs={'Title':'Smoothed (200 segment) along track sea ice thickness estimates from ICESat-2 and ancillary data including ATL10 freeboards and NESOSIM snow loading.',
	'Author':'Alek Petty, Nathan Kurtz, Ron Kwok, Tom Neumann, Thorsten Markus',
	'Contact':'Alek Petty (alek.a.petty@nasa.gov)',
	'Reference':'Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann (2020). Winter Arctic sea ice thickness from ICESat-2 freeboards. J. Geophys. Res. Oceans.'}
	IS2data_short.to_netcdf(outPath+outStr+'_along_track_short.nc')
	
	#---- Plot Data-----
	plot_along_track(IS2data, IS2data_short, outStr)





if __name__ == '__main__':
	


	beam='bnum1'
	releaseStr='rel002'
	runStr='run12'
	cols=['delta_time', 'along_track_distance', 'lat', 'lon', 'ssh_flag', 'seg_length', 'freeboard','snow_depth_N', 'snow_depth_NPdist', 'snow_depth_W99mod5', 'snow_depth_W99mod7', 'snow_depth_W99mod5r', 'snow_depth_W99mod5dist', 'snow_depth_W99mod7dist', 'snow_depth_W99mod5rdist', 'snow_density_N', 'snow_density_W99', 'snow_density_W99r', 'ice_thickness_NPdist', 'ice_thickness_W99mod5dist', 'ice_thickness_unc', 'ice_thickness_uncsys', 'ice_thickness_uncrandom', 'ice_type']

	global outPath
	
	baseDataPath='/cooler/scratch1/aapetty/DataOutput/IS2/'
	dataPath=baseDataPath+'/'+releaseStr+'/'+runStr+'/raw/'
	outPath=baseDataPath+'/'+releaseStr+'/'+runStr+'/along_track/'	

	if not os.path.exists(outPath):
		os.makedirs(outPath)

	#The -01 is for Northern Hemisphere datafiles (-02 is for Southern Hemisphere)

	#thicknessfiles0 = glob(dataPath+'/IS2ATL10-01_2018'+'*'+beam+'*.nc')
	#thicknessfiles1 = glob(dataPath+'/IS2ATL10-01_201901'+'*'+beam+'*.nc')
	#thicknessfiles2 = glob(dataPath+'/IS2ATL10-01_201902'+'*'+beam+'*.nc')
	#thicknessfiles3 = glob(dataPath+'/IS2ATL10-01_201903'+'*'+beam+'*.nc')
	thicknessfiles = glob(dataPath+'/IS2ATL10-01_201904'+'*'+beam+'*.nc')

	#thicknessfiles=thicknessfiles0+thicknessfiles1+thicknessfiles2+thicknessfiles3+thicknessfiles4
	#print('ATL10 files:', ATL10path+'/ATL10-01_'+dateStr+'*.h5')
	print('ATL10 release:',releaseStr)
	print('Processing run:',runStr)
	print('Number of ATL10 files: '+str(np.size(thicknessfiles)))
	#print(thicknessfiles)

	# If you don't want to use concurrent futures and just run over one cpu then use this code
	for thicknessfile in thicknessfiles:
		main(thicknessfile)
	
	#with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:

		# args=((campaign, beam) for beam in beams)
		# print(args)
		# itertools.repeat to add a fixed argument
		# Not very elegant but whatever..d
		#result1=executor.map(main, thicknessfiles)



