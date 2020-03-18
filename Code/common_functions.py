""" common_functions.py
    
    Common functions used in the ICESat-2 sea ice thickness scripts.
    Initial code written by Alek Petty (01/06/2019)
    
    Python dependencies:
        See below for the relevant module imports. Of note:
        xarray/pandas
        netCDF4
        matplotlib
        basemap
        common_functions function library
        Optional (for naive parralel processing):
            itertools
            concurrent.futures
    More information on installation is given in the README file.

    Update history:
        01/06/2020: Version 1.
    
"""


# Load Python modules
import matplotlib
from mpl_toolkits.basemap import Basemap, shiftgrid
from pylab import *
import numpy.ma as ma
from glob import glob
import math
import numpy as np
import xarray as xr
import pandas as pd
import h5py
from scipy.interpolate import griddata
from netCDF4 import Dataset
import dask.array as da
from scipy import stats
from numba import jit, prange


def getIS2gridded(savePathT, outStringT, mapProj, variable='ice_thickness', poleHole=88):
    """ Read in gridded ICESat-2 data
    
    Args:
        savePathT (str): gridded data path
        outStringT (str): label
        mapProj (basemap): map projection

    """

    print(savePathT+'IS2ATL10_'+outStringT+'.nc')
    dIS2 = xr.open_dataset(savePathT+'IS2ATL10_'+outStringT+'.nc')

    varIS2 = array(dIS2[variable])
    latsIS2 = array(dIS2.latitude)
    lonsIS2 = array(dIS2.longitude)

    varIS2 =ma.masked_where(latsIS2>poleHole, varIS2)

    xptsT, yptsT = mapProj(lonsIS2, latsIS2)
    return xptsT, yptsT, lonsIS2, latsIS2, varIS2 

def getIS1gridded(savePathT, outStringT, mapProj,poleHole=85.5):
    """ Read in gridded ICESat-2 data
    
    Args:
        savePathT (str): gridded data path
        outStringT (str): label
        mapProj (basemap): map projection

    """

    print(savePathT+'IS1_'+outStringT+'.nc')
    dIS1 = xr.open_dataset(savePathT+'IS1_'+outStringT+'.nc')

    thicknessIS1 = array(dIS1.ice_thickness)

    latsIS1 = array(dIS1.latitude)
    lonsIS1 = array(dIS1.longitude)
    thicknessIS1 =ma.masked_where(latsIS1>poleHole, thicknessIS1)

    xptsT, yptsT = mapProj(lonsIS1, latsIS1)
    return xptsT, yptsT, lonsIS1, latsIS1, thicknessIS1

def getCS2gsfc(yearStr, mStr):
    """ Read in GSFC CryoSat-2 sea ice thickness data
    
    Args:
        yearStr (str): year
        mStr (str): month

    """

    f = Dataset(dataPathCS2+'/'+yearStr+'/RDEFT4_'+yearStr+mStr+'15.nc', 'r')
    thicknessCS = f.variables['sea_ice_thickness'][:]
    thicknessCS=ma.masked_where(thicknessCS<0, thicknessCS)
    thicknessCS=ma.masked_where(np.isnan(thicknessCS), thicknessCS)

    latsCS = f.variables['lat'][:]
    lonsCS = f.variables['lon'][:]

    xptsT, yptsT = mapProj(lonsCS, latsCS)

    return xptsT, yptsT, thicknessCS
    
    
def getSnowandConverttoThickness(dF, snowDepthVar='snowDepth', snowDensityVar='snowDensity', outVar='iceThickness', rhoi=3):
    """ Grid using nearest neighbour the NESOSIM snow depths to the high-res ICESat-1 freeboard locations"""
    
    # Convert freeboard to thickness
    # Need to copy arrays or it will overwrite the pandas column!
    freeboardT=np.copy(dF['freeboard'].values)
    snowDepthT=np.copy(dF[snowDepthVar].values)
    snowDensityT=np.copy(dF[snowDensityVar].values)
    

    iceDensityT=np.ones_like(dF['freeboard'].values)*916.
    
    if 'ice_density_1' not in dF:
        dF['ice_density_1'] = pd.Series(iceDensityT, index=dF.index)

    if (rhoi==2):
        
        # set lower density for MYI
        iceType=np.copy(dF['ice_type'].values)
        iceDensityT[where(iceType==1)]=882.
        if 'ice_density_2' not in dF:
            dF['ice_density_2'] = pd.Series(iceDensityT, index=dF.index)

    elif (rhoi==3):
        # set lower density for MYI

        iceType=np.copy(dF['ice_type'].values)
        iceDensityT[where(iceType==1)]=899.
        if 'ice_density_3' not in dF:
            dF['ice_density_3'] = pd.Series(iceDensityT, index=dF.index)

    ice_thickness = freeboard_to_thickness(freeboardT, snowDepthT, snowDensityT, iceDensityT)
    #print(ice_thickness)
    dF[outVar] = pd.Series(array(ice_thickness), index=dF.index)
   
    return dF

def getThicknessUncertainty(dF, snowDepthVar, snowDensityVar, iceDensityVar, outVar):
    """ Grid using nearest neighbour the NESOSIM snow depths to the high-res ICESat-1 freeboard locations"""
    
    # Convert freeboard to thickness
    # Need to copy arrays or it will overwrite the pandas column!
    freeboard=np.copy(dF['freeboard'].values)

    
    snow_depth=np.copy(dF[snowDepthVar].values)
    snow_density=np.copy(dF[snowDensityVar].values)
    ice_density=np.copy(dF[iceDensityVar].values)
    
    # Define density values
    water_density=1024.
    #ice_density=915.

    # Add 0.02 to represent precision estimated over flat surfaces
    freeboard_unc=np.copy(dF['freeboard_sigma'].values)+0.02
    #freeboard_unc=0.1*freeboard
    snow_depth_unc=0.05
    snow_density_unc=40

    water_density_unc=0.5
    ice_density_unc=10.

    freeboard_thickness_unc = (freeboard_unc**2)* \
        (water_density/(water_density-ice_density))**2
    
    snow_depth_thickness_unc = (snow_depth_unc**2)* \
        ((snow_density-water_density)/(water_density-ice_density))**2
    
    snow_density_thickness_unc = (snow_density_unc**2)* \
        (snow_depth/(water_density-ice_density))**2
    
    water_density_thickness_unc = (water_density_unc**2)* \
        (((freeboard-snow_depth)/(water_density-ice_density))+
                                    (((snow_depth*water_density)-(freeboard*water_density)-(snow_depth*snow_density))/(water_density-ice_density)**2))**2
    
    ice_density_thickness_unc = (ice_density_unc**2)* \
        (((freeboard*water_density)+(snow_depth*snow_density)-(snow_depth*water_density))/(water_density-ice_density)**2)**2

    # Dropped water density uncertainity as negligible
    rand_uncertainty_squared = (freeboard_thickness_unc+ \
                            snow_depth_thickness_unc + \
                            snow_density_thickness_unc + \
                            ice_density_thickness_unc)

    rand_uncertainty=array(sqrt(rand_uncertainty_squared))

    snow_depth_unc_sys=dF[['snow_depth_NPdist', 'snow_depth_W99mod5rdist', 'snow_depth_Kdist', 'snow_depth_W99mod5dist', 'snow_depth_W99mod7dist']].std(axis=1)
    print(dF.head())
    print('snow depth unc', snow_depth_unc_sys)
    snow_depth_thickness_unc_sys = (snow_depth_unc_sys**2)* \
        ((snow_density-water_density)/(water_density-ice_density))**2

    snow_density_unc_sys=dF[['snow_density_W99','snow_density_W99r', 'snow_density_N']].std(axis=1)
    snow_density_thickness_unc_sys = (snow_density_unc_sys**2)* \
        (snow_depth/(water_density-ice_density))**2

    ice_density_unc_sys=dF[['ice_density_1','ice_density_2', 'ice_density_3']].std(axis=1)
    ice_density_thickness_unc_sys = (ice_density_unc_sys**2)* \
        (((freeboard*water_density)+(snow_depth*snow_density)-(snow_depth*water_density))/(water_density-ice_density)**2)**2


    sys_uncertainty_squared = (snow_depth_thickness_unc_sys + \
                            snow_density_thickness_unc_sys + \
                            ice_density_thickness_unc_sys)
    
    sys_uncertainty=array(sqrt(sys_uncertainty_squared))

    dF['freeboard_unc'] = pd.Series(freeboard_unc, index=dF.index)
    dF[outVar+'random'] = pd.Series(rand_uncertainty, index=dF.index)
    dF[outVar+'sys'] = pd.Series(sys_uncertainty, index=dF.index)
    dF[outVar] = pd.Series(sqrt((rand_uncertainty**2)+(sys_uncertainty**2)), index=dF.index)
   
    return dF


def getATL10FreeboardShotData(fileT, mapProj, beamNum=1, hem='nh'):
    """
    Load ATL10 freeboard shot data from a given beam
    # Store in a Pandas dataframe

    Args:
        freeboardFileT (file): file path of ICESat freeboard data
        mapProj (basemap instance): basemap map projection
        beamNum (int): beam number (1 3, 5 are always strong, the string changes base don IS2 orientation)
    
    Returns:
        dF (var): Dataframe containing freeboard, year, month, day, lon, lat, x, y

    I think the dates are indexed starting from 1 - i.e. month of 1 = January 
        
    """

    print('ATL10 file:', fileT)
    try:
        f1 = h5py.File(fileT, 'r')
        print('got file')
    except:
        return None

    orientation_flag=f1['orbit_info']['sc_orient'][0]
    print('orientation flag:', orientation_flag)
    if (orientation_flag==0):
        print('Backward orientation')
        beamStrs=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
                
    elif (orientation_flag==1):
        print('Forward orientation')
        beamStrs=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
        
    elif (orientation_flag==2):
        print('Transitioning, do not use for science!')
    
    beamStr=beamStrs[beamNum-1]
    print(beamStr)

    #freeboard_confidence=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['beam_fb_confidence'][:]
    #freeboard_quality=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['beam_fb_quality_flag'][:]
    try:
        freeboard=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['beam_fb_height'][:]
    except:
        print ('no valid freeboard data for that beam')
        return None

    segment_length=f1[beamStr]['freeboard_beam_segment']['height_segments']['height_segment_length_seg'][:]
      
    freeboard_sigma=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['beam_fb_sigma'][:]
    freeboard_sigma[where(freeboard_sigma>0.1)]=0.1
    
    ssh_flag=f1[beamStr]['freeboard_beam_segment']['height_segments']['height_segment_ssh_flag'][:]

    lons=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['longitude'][:]
    lats=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['latitude'][:]
    deltaTime=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['delta_time'][:]
    algTrkDist=f1[beamStr]['freeboard_beam_segment']['beam_freeboard']['seg_dist_x'][:]

    dF = pd.DataFrame({'freeboard':freeboard, 'freeboard_sigma':freeboard_sigma, 'seg_length':segment_length, 'ssh_flag':ssh_flag, 'along_track_distance':algTrkDist, 'lon':lons, 'lat':lats, 'delta_time':deltaTime})

    if (hem=='nh'):
        print('Keeping Northern lat values only')
        dF = dF[(dF['lat']>0)]
    else:
        dF = dF[(dF['lat']<0)]
    #dF = dF[(dF['freeboard']>0)]
    dF = dF[(dF['freeboard']<10)]

    # Decide here if we want to also filter based on the confidence and/or quality flag

    # Reset row indexing
    dF=dF.reset_index(drop=True)
    #print(dF.head(3))
    
    # Change this to read in from file not the string
    fileStr=fileT.split("/ATL10-")[-1]
    fileStr=fileStr.split("_")[1]
    print('fileStr', fileStr)
    dF['year'] = int(fileStr[0:4])
    dF['month'] = int(fileStr[4:6])
    dF['day'] = int(fileStr[6:8])

    print('Year:', fileStr[0:4], 'Month:', fileStr[4:6], 'Day:', fileStr[6:8])
    

    xpts, ypts=mapProj(dF['lon'].values, dF['lat'].values)

    dF['xpts'] = pd.Series(xpts, index=dF.index)
    dF['ypts'] = pd.Series(ypts, index=dF.index)

    dF['beamStr'] = beamStr
    dF['beamNum'] = str(beamNum)
    
    return dF, fileStr


def getProcessedATL10Shotdata(dataOutPathM, runStrT, campaignStrT, cols, yearStr='2018', monStr='*', dayStr='*', fNum=-1, beamStr='gt1r'):
    """
    Load ICESat-2 thickness data produced from the raw ATL10 shot data
        
    """
    files=glob(dataOutPathM+runStrT+'/'+campaignStrT+'/IS2ATL10*'+yearStr+monStr+dayStr+'*'+beamStr)
    print(dataOutPathM+runStrT+'/'+campaignStrT+'/IS2ATL10*'+yearStr+monStr+dayStr+'*'+beamStr)
    print('Number of files:', size(files))

    if (fNum>=0):
        # if we want to just pick one of the files available that day...
        IS2dataAll = pd.read_pickle(files[fNum])
        if (cols!='all'):
            IS2dataAll=IS2dataAll[cols]
        return IS2dataAll


    if size(files)>1:
        IS2dataAll = pd.read_pickle(files[0])
        IS2dataAll=IS2dataAll[cols]
        for x in range(1, size(files)):
            IS2dataT = pd.read_pickle(files[x])
            if (cols!='all'):
                IS2dataT=IS2dataT[cols]
            IS2dataAll = pd.concat([IS2dataT, IS2dataAll])
        
            print ('File:', str(x+1)+'/'+str(size(files)), IS2dataAll.shape)
        IS2dataAll=IS2dataAll.reset_index(drop=True)
        
        return IS2dataAll
    elif size(files)==1:
        IS2dataAll = pd.read_pickle(files[0])
        if (cols!='all'):
            IS2dataAll=IS2dataAll[cols]
        return IS2dataAll
    else:
        print('No files')
        return



def getProcessedATL10ShotdataNCDF(dataPathT, yearStr='*', monStr='*', dayStr='*', fNum=-1, concat=0, ssh_mask=0, minseg=0, maxseg=0, beamStr='gt1r', vars=[], smoothingWindow=0):
    """
    Load ICESat-2 thickness data produced from the raw ATL10 shot data
        
    """
    
    print(dataPathT+'IS2ATL10-01_'+yearStr+monStr+dayStr+'*'+'_*'+beamStr+'*.nc')
    files=glob(dataPathT+'IS2ATL10-01_'+yearStr+monStr+dayStr+'*'+'_*'+beamStr+'*.nc')
    print('Number of files:', size(files))
    #import pdb; pdb.set_trace()
    print(ssh_mask, minseg, maxseg)
    #testFile = Dataset(files[0])
    #print(testFile.variables.keys())
    if (fNum>-0.5):
        if (size(vars)>0):
            IS2dataAll= xr.open_dataset(files[fNum], engine='h5netcdf')[vars]
            print(IS2dataAll)
        else:
            
            IS2dataAll= xr.open_dataset(files[fNum], engine='h5netcdf')

    else:

        if (concat==1):
            IS2files=[xr.open_dataset(file,chunks={}, engine='h5netcdf')[vars] for file in files] 
            darray = xr.concat(IS2files, dim='index')
            return darray

        # apparently autoclose assumed so no longer need to include the True flag
        if (size(vars)==1):
            IS2files=[xr.open_dataset(file,chunks={}, engine='h5netcdf')[vars] for file in files] 
            darray = xr.concat(IS2files, dim='index')
            #darray=darray.chunk(f3.size/40)
            return darray
        if (size(vars)>1):
            print('variable')
            #files=glob(dataPathT+'/IS2ATL10*'+yearStr+monStr+dayStr+'*'+'_'+beamStr+'*.nc')
            #IS2files=[xr.open_dataset(file, engine='h5netcdf')[vars] for file in files] 
            #IS2dataAll=xr.concat(IS2files, 'index')
            IS2dataAll= xr.open_mfdataset(dataPathT+'/IS2ATL10-01_'+yearStr+monStr+dayStr+'*'+'_'+beamStr+'*.nc', concat_dim='index', parallel=True)[vars]
        else:
            IS2dataAll= xr.open_mfdataset(dataPathT+'/IS2ATL10-01_'+yearStr+monStr+dayStr+'*'+'_'+beamStr+'*.nc', parallel=True)
        #IS2dataAll = pd.read_pickle(files[0])
    print(IS2dataAll.info)
    
    #IS2dataAll=IS2dataAll[vars]
    #print(IS2dataAll.info)

    if (ssh_mask==1):
        print('removing data where ssh_flag = 1')
        IS2dataAll=IS2dataAll.where(IS2dataAll.ssh_flag<0.5, drop=True)
        print(IS2dataAll.info)

    if (minseg>0):
        print('seg range:', np.amin(IS2dataAll['seg_length'].values), np.max(IS2dataAll['seg_length'].values))
        print('min seg filter..', minseg)
        IS2dataAll=IS2dataAll.where(IS2dataAll['seg_length']>minseg, drop=True)
        print('seg range:', np.amin(IS2dataAll['seg_length'].values), np.max(IS2dataAll['seg_length'].values))
    
    if (maxseg>0):
        print('max seg filter..', maxseg)
        IS2dataAll=IS2dataAll.where(IS2dataAll['seg_length']<maxseg, drop=True)
        print('seg range:', np.amin(IS2dataAll['seg_length'].values), np.max(IS2dataAll['seg_length'].values))

    if (smoothingWindow>0.5):
        print('smoothing...')
        # If we want to smooth the datasets
        seg_length=IS2dataAll['seg_length']
        
        seg_weightedvarR=seg_length.rolling(index=smoothingWindow, center=True).mean()
        seg_weightedvar=seg_weightedvarR[int(smoothingWindow/2):-int(smoothingWindow/2):smoothingWindow]
       # print(seg_weightedvar)

        ds = seg_weightedvar.to_dataset(name = 'seg_length')
        #seg_weightedvars.append(seg_weightedvar)
        # Skip the first one as that's always (should be) the seg_length
        for var in vars[1:]:
            print('Coarsening'+var+'...')
            varIS2=IS2dataAll[var]
            seg_weightedvarR=varIS2*seg_length.rolling(index=smoothingWindow, center=True).sum()/seg_length.rolling(index=smoothingWindow, center=True).sum()
            seg_weightedvar=seg_weightedvarR[int(smoothingWindow/2):-int(smoothingWindow/2):smoothingWindow] 
            #print(seg_weightedvar)
            ds[var] = seg_weightedvar

            #seg_weightedvars.append(seg_weightedvar)
            print('Coarsened var')
        #Merge the coarsened arrays
        #seg_weightedvarsM=xr.merge(seg_weightedvars)
        ds=ds.reset_index('index', drop=True)

        #print('Rechunking...')
        #ds=ds.chunk(2000)
        #print('Rechunked')
        print(ds)

        IS2dataAll.close()
        return ds
    else:

        return IS2dataAll

def getProcessedIS1(dataPathT, campaignStr, vars=[], fNum=-1, smoothingWindow=0):
    """
    Load ICESat-2 thickness data produced from the raw ATL10 shot data
        
    """
    
    print(dataPathT+'IS1*'+campaignStr+'.nc')
    files=glob(dataPathT+'IS1*'+campaignStr+'.nc')
    print('Number of files:', size(files))
    #import pdb; pdb.set_trace()

    #testFile = Dataset(files[0])
    #print(testFile.variables.keys())
    if (fNum>-0.5):
        if (size(vars)>0):
            IS1dataAll= xr.open_dataset(files[fNum], engine='h5netcdf')[vars]
            print(IS1dataAll)
        else:
            
            IS1dataAll= xr.open_dataset(files[fNum], engine='h5netcdf')

    else:
        # apparently autoclose assumed so no longer need to include the True flag
        if (size(vars)>0):
            IS1dataAll= xr.open_mfdataset(dataPathT+'/IS1*'+campaignStr+'.nc', engine='h5netcdf', data_vars=vars, parallel=True)
        else:
            IS1dataAll= xr.open_mfdataset(dataPathT+'/IS1*'+campaignStr+'.nc', engine='h5netcdf', parallel=True)
        #IS2dataAll = pd.read_pickle(files[0])
    print(IS1dataAll.info)
    


    if (smoothingWindow>0):
       
        
        for var in vars:
            print('Coarsening'+var+'...')
            varIS2=IS1dataAll[var]
            svarR=varIS1.rolling(index=smoothingWindow, center=True).mean()
            if (var==vars[0]):
                ds = svarR.to_dataset(name = var)
            else:
                ds[var] = svarR

            print('Coarsened var')
        #Merge the coarsened arrays
        ds=ds.reset_index('index', drop=True)

        print(ds)
        return ds
    else:

        return IS1dataAll


def reset_matplotlib():
    """
        Reset matplotlib to a common default.
    
    """
    
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    rcParams['ytick.major.size'] = 2
    rcParams['axes.linewidth'] = .6
    rcParams['lines.linewidth'] = .6
    rcParams['patch.linewidth'] = .6
    rcParams['ytick.labelsize']=8
    rcParams['xtick.labelsize']=8
    rcParams['legend.fontsize']=9
    rcParams['font.size']=9
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def getNesosimDates(dF, snowPathT):
    """ Get dates from NESOSIM files"""

    # This will come from the distinct rows of the IS-1 data eventually, but for now the data only span a day or two, so not much change in snow depth..
    dayS=dF['day'].iloc[0]
    monthS=dF['month'].iloc[0]
    monthF=dF['month'].iloc[-1]
    yearS=dF['year'].iloc[0]
    dateStr= getDate(dF['year'].iloc[0], dF['month'].iloc[0], dF['day'].iloc[0])

    print ('Date:', yearS, monthS, dayS)
    # Find the right NESOSIM data file based on the freeboard dates
    if (monthS<6):
        fileNESOSIM = glob(snowPathT+'*'+str(yearS)+'.nc')[0]
    else:
        fileNESOSIM = glob(snowPathT+'*'+str(yearS)+'*.nc')[0]

    if (monthS>5 & monthF==5):
        print ('WARNING! LACK OF SNOW DATA')

    return fileNESOSIM, dateStr

def plotMap4(dF, mapProj, figPathT, outStr, vars=['freeboard', 'snowDepthN', 'snowDensityN', 'iceThicknessN']):
    """ plot data on 4 maps as a check"""

    rcParams['axes.labelsize'] = 8
    rcParams['xtick.labelsize'] = 8
    rcParams['ytick.labelsize'] = 8
    rcParams['font.size'] = 8
    rcParams['lines.linewidth'] = 1
    rcParams['patch.linewidth'] = 1
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    cbarx=[0.03, 0.56, 0.03, 0.56]
    cbary=[0.55, 0.55, 0.07, 0.07]
    vmins=[0, 0, 150, 0]
    vmaxs=[1.5, 0.4, 380, 5]
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(6, 7))

    for x in range(4):
        ax=axs.flatten()[x]
        sca(ax)

        im1=hexbin(dF['xpts'].values, dF['ypts'].values, C=dF[vars[x]].values, gridsize=100, 
               cmap=cm.viridis, vmin=vmins[x], vmax=vmaxs[x], zorder=2, rasterized=True)
        #m.fillcontinents(color='',lake_color='grey', zorder=3)
        mapProj.drawcoastlines(linewidth=0.25, zorder=5)
        mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
        mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)

        cax = fig.add_axes([cbarx[x], cbary[x], 0.4, 0.025])

        #cax = fig.add_axes([0.25, 0.15, 0.5, 0.04])
        cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
        cbar.set_label(vars[x]+' ', labelpad=3, fontsize=14)

    subplots_adjust(left = 0.01, right = 0.99, bottom=0.08, top = 0.99, hspace=0.1)
    plt.savefig(figPathT+'/test'+outStr+'4.png')
    plt.close('all')

def correlateVars(var1, var2):
    """correlate two variables"""

    trend, intercept, r_a, prob, stderr = stats.linregress(var1, var2)
    sig = 100.*(1.-prob)
    return trend, sig, r_a, intercept 


def getOSISAFconcDay(osisafPath, mProj, yearT, month, day, mask=1, maxConc=0, lowerConc=1, fillPhole=1, regionMask=1):
    """
    Get daily ice conc from OSISAF

    month index starts at 1
    day index starts at 1

    """
    pmask=get_pmask(yearT, month)

    mStr='%02d' %(month)
    dStr='%02d' %(day) 
    fileT=osisafPath+str(yearT)+'/'+mStr+'/ice_conc_nh_polstere-100_multi_'+str(yearT)+mStr+dStr+'1200.nc'
            

    f = Dataset(fileT, 'r')
    print(fileT)

    # read lat/lon (start points) and lat1/lon1 (end points)
    lon = (f.variables['lon'][:])
    lat = (f.variables['lat'][:])

    # transform to map project coordinates (basemap's axes, assume they have unit of m)
    x0, y0=mProj(lon, lat)

    # Convert percent to conc!
    conc = (f.variables['ice_conc'][0])/100.
    
    f.close()

    if (mask==1):
        conc = ma.masked_where(conc>1., conc)
    
    if (maxConc==1):
        conc = ma.where(conc>1.,0, conc)

    if (lowerConc==1):
        conc = ma.where(conc<0.15,0, conc)

    if (fillPhole==1):
        concHole=ma.mean(conc[(lat>pmask-2.) & (lat<pmask-1.5)])
        #print concHole
        conc = where((lat >=pmask-2.), concHole, conc)

    if (regionMask==1):
        region_mask, xptsI, yptsI = get_region_mask_sect('../../AncData/', mProj, xypts_return=1)

        conc[where(region_mask>18)]=0

    #print iceConcDay
    #iceConcDay=ma.filled(np.nan)
    
    return conc, x0, y0

def turningpoints(yvals, xvals):
    dx = np.diff(yvals)
    dx2=np.diff(dx)
    tpts=where((dx[1:] * dx[:-1] < 0)&(dx2<0))
    tpts=[t+1 for t in tpts][0]
    print(tpts)
    print(yvals[tpts])
    # limit these to values above a threshold to avoid weird bins at the tail
    #idx=yvals[tpts]>0.1*np.amax(yvals)
    modes=[]
    for tpt in tpts:
        print(yvals[tpt])
        if (yvals[tpt]>0.1*np.amax(yvals)):
            modes.append(xvals[tpt])

    return dx, dx2, modes

def get_region_mask_sect_labels(region):

    labels=['non-region oceans', 'Sea of Okhotsk and Japan', 'Bering Sea', 'Hudson Bay','Gulf of St. Lawrence', \
    'Baffin Bay/D. Strait/Lab. Sea', 'Greenland Sea', 'Barents Sea', 'Kara Sea','Laptev Sea', 'E. Siberian Sea',\
    'Chukchi Sea', 'Beaufort Sea', 'Canadian Archipelago', 'Central Arctic', 'Land', 'Coast']
    
    return labels[region-1]

def get_cdr_conc(concdataPath, mplot, yearStr, monStr):
    fileT=glob(concdataPath+'*'+yearStr+monStr+'*.nc')[0]
    fconc = Dataset(fileT, 'r')
    iceConc = fconc.variables['seaice_conc_monthly_cdr'][0]
    clat = fconc.variables['latitude'][:]
    clon = fconc.variables['longitude'][:]
    xptsc, yptsc=mplot(clon, clat)
    return xptsc, yptsc, iceConc

def get_region_mask_sect(datapath, mplot, xypts_return=0):
    datatype='uint8'
    file_mask = datapath+'/sect_fixed_n.msk'
    # 1   non-region oceans
    # 2   Sea of Okhotsk and Japan
    # 3   Bering Sea
    # 4   Hudson Bay
    # 5   Gulf of St. Lawrence
    # 6   Baffin Bay/Davis Strait/Labrador Sea
    # 7   Greenland Sea
    # 8   Barents Seas
    # 9   Kara
    # 10   Laptev
    # 11   E. Siberian
    # 12   Chukchi
    # 13   Beaufort
    # 14   Canadian Archipelago
    # 15   Arctic Ocean
    # 20   Land
    # 21   Coast
    fd = open(file_mask, 'rb')
    region_mask = fromfile(file=fd, dtype=datatype)
    region_mask = reshape(region_mask, [448, 304])

    #xpts, ypts = mplot(lons_mask, lats_mask)
    if (xypts_return==1):
        mask_latf = open(datapath+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(datapath+'/psn25lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask

def get_pmask(year, month):
    # Get the pole hole latitude for the various passive microwave sensors
    # Remove half a degree as gridding around the pole hole edge
    if (year<1987):
        pmask=84.4
    elif((year==1987)&(month<=6)):
        pmask=84.4
    elif ((year==1987)&(month>6)):
        pmask=86.7
    elif ((year>1987)&(year<2008)):
        pmask=87.2
    else:
        pmask=89.2
    
    return pmask

def get_psnlatslons(data_path, res=25):
    if (res==25):
        # 25 km grid
        mask_latf = open(data_path+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn25lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])
    elif (res==12):
        # 12.5 km grid
        mask_latf = open(data_path+'/psn12lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn12lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [896, 608])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [896, 608])
    elif (res==6):
        # 12.5 km grid
        mask_latf = open(data_path+'/psn06lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn06lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [1792, 1216])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [1792, 1216])

    return lats_mask, lons_mask

def distributeSnow(dF, inputSnowDepth='snowDepthNP', outSnowVar='snowDepthNPdist', consIterations=11, gridSize=100000, version='V3'):
    """
    Do the distribution from the snow data already loaded in

    Args:
        dF (data frame): Pandas dataframe
        mapProj (basemap instance): Basemap map projection
        fileSnow (string): NESOSIM file path
        dateStr (string): date string
        outSnowVar (string): Name of snow depth column
        consIterations (var, default of 11): number of terative loops to try and conserve snow after the freeboard re-distribution
        gridSize (var, default of 100000): along track distance in meters
    Returns:
        dF (data frame): Pandas dataframe updated to include colocated NESOSIM (and dsitributed) snow data

    """
    # num of points required to do a snow redistribution calculation
    numPtsReq=100

    #dN = xr.open_dataset(fileSnow)

    # Get NESOSIM snow depth and density data for that date
    # Could move this into the loop if there is a significant date spread in the freeboard data (appears not to be the case though!).
    # Improves processing speed. 
    
    #dNday = dN.sel(day=int(dateStr))
    # Get NESOSIM coordinates (constant with time)
    #lonsN = array(dNday.longitude)
    #latsN = array(dNday.latitude)
    #xptsN, yptsN = mapProj(lonsN, latsN)

    #dateStr= getDate(dF['year'].iloc[0], dF['month'].iloc[0], dF['day'].iloc[0])
    #print (dateStr)
    #print ('Freeboard (m):', dF['freeboard'].iloc[x])
    #dNday = dN.sel(day=int(dateStr))
    #snowDepthNDay = array(dNday.snowDepth)
    #snowDensityNDay = array(dNday.density)
    #iceConcNDay = array(dNday.iceConc)

    #mask=where((snowDepthNDay>0.01)&(snowDepthNDay<1)&(iceConcNDay>0.01)&np.isfinite(snowDensityNDay))

    #snowDepthNDay = snowDepthNDay[mask]
    #snowDensityNDay = snowDensityNDay[mask]
    #xptsNDay = xptsN[mask]
    #yptsNDay = yptsN[mask]
    xptsT=dF['xpts'].values
    yptsT=dF['ypts'].values
    freeboardsT=dF['freeboard'].values
    
    snowDepthUndistributed=dF[inputSnowDepth].values

    #snowDepthGISs=ma.masked_all(size(freeboardsT))
    #snowDensityGISs=ma.masked_all(size(freeboardsT))
    #for x in range(size(freeboardsT)):
    #   snowDensityGIS = griddata((xptsNDay, yptsNDay), snowDensityNDay, (xptsT[x], yptsT[x]), method='nearest') 
    #    snowDensityGISs[x]=snowDensityGIS

    # consered/distributed snow depth
    snowDepthDistsC=ma.masked_all(size(freeboardsT))
    # Unconsered/distributed snow depth
    #snowDepthDistsUC=ma.masked_all(size(freeboardsT))

    alongTrackDists = sqrt((xptsT-xptsT[0])**2+(yptsT-yptsT[0])**2)
    numGrids=int(ceil(alongTrackDists[-1]/gridSize))
    print('Number of 100km grids in file:', numGrids)
    
    for gridIdx in range(numGrids):
        #print(gridIdx)
        gridPts=where((alongTrackDists>=(gridSize*gridIdx))&(alongTrackDists<(gridSize*(gridIdx+1))))[0]
        #print('grididx', gridIdx, ', gridPts:', gridPts, ', size:', size(gridPts), ', dists:', alongTrackDists[gridPts], int((size(gridPts) - 1)/2))
        
        if (size(gridPts)>100):
            #midIdx = int((size(gridPts) - 1)/2)
            #midGridIdx=gridPts[midIdx]
            #print('grid points:', gridPts, midIdx)
            
            if (version=='V3'):
                snowDepthDistsC[gridPts] = snowDistributionV3(ma.mean(snowDepthUndistributed[gridPts]), freeboardsT[gridPts], ma.mean(freeboardsT[gridPts]), numIter=consIterations)
            elif (version=='V4'):
                snowDepthDistsC[gridPts] = snowDistributionV4(ma.mean(snowDepthUndistributed[gridPts]), freeboardsT[gridPts], ma.mean(freeboardsT[gridPts]), numIter=consIterations)
            #print (snowDepthUndistributed[gridPts])
            #print (snowDepthDistsC[gridPts])
            
            #snowDepthDistsUC[gridPts] = snowDistributionV3(ma.mean(snowDepthUndistributed[gridPts]), freeboardsT[gridPts], ma.mean(freeboardsT[gridPts]), numIter=1)
            
            #print (snowDepthDistsUC[gridPts])

        else:
            # Just find the nearest snow depth
            for gridPt in gridPts:
                snowDepthDistsC[gridPt]=snowDepthUndistributed[gridPt] 
                #snowDepthDistsUC[gridPt]=snowDepthUndistributed[gridPt] 
                
    dF[outSnowVar] = pd.Series(snowDepthDistsC, index=dF.index)
    #dF[outSnowVar+'UC'] = pd.Series(snowDepthDistsUC, index=dF.index)
    return dF


def distributeSnowKwok(dF, inputSnowDepth='snowDepthN', outSnowVar='snowDepthKdist'):
    """
    Do the distribution from the snow data already loaded in

    Args:
        dF (data frame): Pandas dataframe

    Returns:
        dF (data frame): Pandas dataframe updated to include colocated NESOSIM (and dsitributed) snow data
        
    """

    xptsT=dF['xpts'].values
    yptsT=dF['ypts'].values
    freeboardsT=dF['freeboard'].values
    
    snowDepthUndistributed=dF[inputSnowDepth].values

    snowDepthDistsC=ma.masked_all(size(freeboardsT))

    for x in range(size(freeboardsT)):
        frac=sigmoidFunc(freeboardsT[x]/snowDepthUndistributed[x])
        snowDepthDistsC[x]=frac*snowDepthUndistributed[x]

        if (snowDepthDistsC[x] > freeboardsT[x]): 
            # Dont let snow be greater than the freeboard (this is where it ges harder to conserve snow)
            snowDepthDistsC[x] = freeboardsT[x]

    dF[outSnowVar] = pd.Series(snowDepthDistsC, index=dF.index)
    dF.head(3)

    return dF

def sigmoidFunc(x):
    return 1./(1+exp(-5*(x-0.7)))

#@jit
def gridNESOSIMtoFreeboard(dF, mapProj, fileSnow, dateStr, outSnowVar='snowDepthN', outDensityVar='snowDensityN', returnMap=0):
    """
    Load relevant NESOSIM snow data file and assign to freeboard values

    Args:
        dF (data frame): Pandas dataframe
        mapProj (basemap instance): Basemap map projection
        fileSnow (string): NESOSIM file path
        dateStr (string): date string
        outSnowVar (string): Name of snow depth column
        outDensityVar (string): Name of snow density column

    Returns:
        dF (data frame): dataframe updated to include colocated NESOSIM (and dsitributed) snow data

    """

    dN = xr.open_dataset(fileSnow)

    # Get NESOSIM snow depth and density data for that date
    # Should move this into the loop if there is a significant date cahgne in the freeboard data.
    # Not done this to improve processing speed. 
    dNday = dN.sel(day=int(dateStr))
    
    lonsN = array(dNday.longitude)
    latsN = array(dNday.latitude)
    xptsN, yptsN = mapProj(lonsN, latsN)

    # Get dates at start and end of freeboard file
    dateStrStart= getDate(dF['year'].iloc[0], dF['month'].iloc[0], dF['day'].iloc[0])
    dateStrEnd= getDate(dF['year'].iloc[-1], dF['month'].iloc[-1], dF['day'].iloc[-1])
    print('Check dates (should be within a day):', dateStr, dateStrStart, dateStrEnd)
  
    snowDepthNDay = array(dNday.snowDepth)
    snowDensityNDay = array(dNday.density)
    iceConcNDay = array(dNday.iceConc)
    
    # Remove data where snow depths less than 0 (masked).
    # Might need to chek if I need to apply any other masks here.
    mask=where((snowDepthNDay>0.01)&(snowDepthNDay<1)&(iceConcNDay>0.01)&np.isfinite(snowDensityNDay))

    snowDepthNDay = snowDepthNDay[mask]
    snowDensityNDay = snowDensityNDay[mask]
    xptsNDay = xptsN[mask]
    yptsNDay = yptsN[mask]

    # Load into array, sppeds up later computation and may aid parallelization
    freeboardsT=dF['freeboard'].values
    xptsT=dF['xpts'].values
    yptsT=dF['ypts'].values
    
    # I think it's better to declare array now so memory is allocated before the loop?
    snowDepthGISs=ma.masked_all(size(freeboardsT))
    snowDensityGISs=ma.masked_all(size(freeboardsT))
    #snowDepthDists=ma.masked_all(size(freeboardsT))

    #for x in prange(size(freeboardsT)):
    for x in range(size(freeboardsT)):
        
        # Could embed the NESOSIM dates here
        
        # Use nearest neighbor to find snow depth at IS2 point
        #snowDepthGISs[x] = griddata((xptsDay, yptsDay), snowDepthDay, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        #snowDensityGISs[x] = griddata((xptsDay, yptsDay), densityDay, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest')

        # Think this is the much faster way to find nearest neighbor!
        dist = sqrt((xptsNDay-xptsT[x])**2+(yptsNDay-yptsT[x])**2)
        index_min = np.argmin(dist)
        snowDepthGISs[x]=snowDepthNDay[index_min]
        snowDensityGISs[x]=snowDensityNDay[index_min]
        #print(snowDepthNDay[index_min], densityNDay[index_min])
        
    dF[outSnowVar] = pd.Series(snowDepthGISs, index=dF.index)
    dF[outDensityVar] = pd.Series(snowDensityGISs, index=dF.index)

    # SNOW REDISTRIBUTION
    #for x in range(size(freeboardsT)):
        
        # Find the mean freebaord in this vicinitiy
        # ICESat-1 has a shot every 172 m, so around 600 shots = 100 km
    #    meanFreeboard = ma.mean(freeboardsT[x-300:x+300])
    #    snowDepthDists[x] = snowDistribution(snowDepthGISs[x], freeboardsT[x], meanFreeboard)

    #dF[outSnowVar+'dist'] = pd.Series(snowDepthDists, index=dF.index)

    #print ('Snow depth (m): ', snowDepthGIS)
    #print ('Snow density (kg/m3): ', snowDensityGIS)
    #print ('Snow depth (m): ', snowDepthDists)
    
    if (returnMap==1):
        return dF, xptsN, yptsN, dNday, 
    else:
        return dF


def snowDistributionV3(meanSnowDepthT, freeboardsT, meanFreeboardT, numIter=11):
    """
    Snow distribution function
    Written by Alek Petty (02/01/2019)

    Args:
        meanSnowDepthT (var): mean snow depth (m)
        freeboardT (var): along track freeboard (m)
        meanFreeboardT (var): mean freeboard (m)
        
        numIter (var): number of iterations to minimize the functional fit
        
    Returns:
        snowDepthDists (var): along track snow depth distributed to higher resolution
        
    """

    

    fb_cutoff = (0.71*meanSnowDepthT) + (0.20*meanFreeboardT) + 0.146
    hs_thick_adj = 0.0
    hs_thick_prev_adj = 0.0

    #print('Mean snow:', meanSnowDepthT, 'mean freeboard:', meanFreeboardT, 'fcutoff', fb_cutoff)


    for iter in range(numIter): # Set this to zero to not do an iterative loop to conserve snow
        # This doesn't seem to be doing much
        #print ('iter:', iter)
        hs_thick = (1.03*meanSnowDepthT) + 0.0083 + hs_thick_adj 

        distSnow = np.zeros(len(freeboardsT))
        x=0
        for freeboardT in freeboardsT: #len(freeboardT)): #Apparently freeboardT doesn't start with index 0!
            #print ('freeboardT:', freeboardT)
            hs_thin = hs_thick * freeboardT/fb_cutoff    #linear way
            
            distSnow[x] = hs_thin

            if freeboardT >= fb_cutoff: 
                # For large freeboards where there is no correlation between freeboard and snow depth
                distSnow[x] = hs_thick   

            if distSnow[x] > freeboardT: 
                # Dont let snow be greater than the freeboard (this is where it ges harder to conserve snow)
                distSnow[x] = freeboardT
            x+=1
        #print('distsnow:', distSnow)
        
        # Calculate hs_thick_adj
        hs_thick_adj = (meanSnowDepthT - np.mean(distSnow)) + hs_thick_prev_adj
        if abs(hs_thick_adj-hs_thick_prev_adj) < 0.01:
            break 
        hs_thick_prev_adj = hs_thick_adj
    
    #print('Final dist snow:', distSnow)
    return distSnow


def snowDistributionV4(meanSnowDepthT, freeboardsT, meanFreeboardT, numIter=11):
    """
    Snow distribution function
    Written by Alek Petty (02/01/2019)

    Args:
        meanSnowDepthT (var): mean snow depth (m)
        freeboardT (var): along track freeboard (m)
        meanFreeboardT (var): mean freeboard (m)
        
        numIter (var): number of iterations to minimize the functional fit
        
    Returns:
        snowDepthDists (var): along track snow depth distributed to higher resolution
        
    """
    

    adjSnow = np.zeros(len(freeboardsT))
    x=0
    for freeboardT in freeboardsT: #len(freeboardT)): #Apparently freeboardT doesn't start with index 0!
        # work out where we should have got rid of some snow in the accumulation model because of unphysical (higher than freeboard) values
        if meanSnowDepthT > freeboardT: 
            # Dont let snow be greater than the freeboard (this is where it ges harder to conserve snow)
            adjSnow[x] = freeboardT
        else:
            adjSnow[x]=meanSnowDepthT
        x+=1
    
    adjmeanSnowDepthT=np.mean(adjSnow)

    fb_cutoff = (0.7*meanSnowDepthT) + (0.2*meanFreeboardT) + 0.15
    hs_thick_adj = 0.0
    hs_thick_prev_adj = 0.0

    #print('Mean snow:', meanSnowDepthT, 'Adj mean snow:', adjmeanSnowDepthT, 'mean freeboard:', meanFreeboardT, 'fcutoff', fb_cutoff)


    for iter in range(numIter): # Set this to zero to not do an iterative loop to conserve snow
        # This doesn't seem to be doing much
        #print ('iter:', iter)
        hs_thick = adjmeanSnowDepthT + hs_thick_adj 

        distSnow = np.zeros(len(freeboardsT))
        x=0
        for freeboardT in freeboardsT: #len(freeboardT)): #Apparently freeboardT doesn't start with index 0!
            #print ('freeboardT:', freeboardT)
            hs_thin = hs_thick * freeboardT/fb_cutoff    #linear way
            
            distSnow[x] = hs_thin

            if freeboardT >= fb_cutoff: 
                # For large freeboards where there is no correlation between freeboard and snow depth
                distSnow[x] = hs_thick   

            if distSnow[x] > freeboardT: 
                # Dont let snow be greater than the freeboard (this is where it ges harder to conserve snow)
                distSnow[x] = freeboardT
            x+=1
        #print('distsnow:', distSnow)
        
        # Calculate hs_thick_adj
        hs_thick_adj = (adjmeanSnowDepthT - np.mean(distSnow)) + hs_thick_prev_adj
        if abs(hs_thick_adj-hs_thick_prev_adj) < 0.01:
            break 
        hs_thick_prev_adj = hs_thick_adj
    
    #print('Final dist snow:', distSnow)
    return distSnow


#@jit
def assignRegionMask(dF, mapProj, ancDataPath='../../AncData/'):
    """
    Apply NSIDC sectional region mask to the ICESat-2 data

    # 1   non-region oceans
    # 2   Sea of Okhotsk and Japan
    # 3   Bering Sea
    # 4   Hudson Bay
    # 5   Gulf of St. Lawrence
    # 6   Baffin Bay/Davis Strait/Labrador Sea
    # 7   Greenland Sea
    # 8   Barents Seas
    # 9   Kara
    # 10   Laptev
    # 11   E. Siberian
    # 12   Chukchi
    # 13   Beaufort
    # 14   Canadian Archipelago
    # 15   Arctic Ocean
    # 20   Land
    # 21   Coast

    Args:
        dF (data frame): original data frame
        mapProj (basemap instance): basemap map projection
        
        
    Returns:
        dF (data frame): data frame including ice type column (1 = multiyear ice, 0 = everything else)

    """


    region_mask, xptsI, yptsI = get_region_mask_sect(ancDataPath, mapProj, xypts_return=1)

    xptsI=xptsI.flatten()
    yptsI=yptsI.flatten()
    region_mask=region_mask.flatten()

    #iceTypeGs=[]
    regionFlags=ma.masked_all((size(dF['freeboard'].values)))
    for x in range(size(dF['freeboard'].values)):
        # Find nearest region
        dist=sqrt((xptsI-dF['xpts'].iloc[x])**2+(yptsI-dF['ypts'].iloc[x])**2)
        index_min = np.argmin(dist)
        regionFlags[x]=int(region_mask[index_min])
        # This is waht I normally do but it's much slower!
        # I checked and they gave the same answers
        # iceTypeG2 = griddata((xpts_type, ypts_type), ice_typeT2, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        # print(iceTypeG)
        # iceTypeGs.append(iceTypeG)

    dF['region_flag'] = pd.Series(regionFlags, index=dF.index)

    return dF


def getIceTypeRaw(iceTypePathT, mapProj, dayStr, monStr, yearStr, res=1):
    """
    Get ice type based on the OSI-SAF ice type product, 
    extended back tot he SSMI period as part of the Copernicus Climate Change (C3S) service.
    Thanks to Thomas Lavergne and Singe Aaboe 
    
    Data from 2006 onwards is from OSI-SAF (Lavergene et al.)
    Data prior to 2006 is produced through funding provided by the Copernicus CDR (Aaboe et al.,)
    Product User Guide and Specification (PUGS): Sea Ice Edge and Type version 1, Algorithm Theoretical Basis Document (ATBD): Sea Ice Edge and Type version 2";
    
    Assumes minimal change in ice type across a day so just use the first row for date information and apply to all rows

    Returns:
        x, y, icetype

    """

    dateStr=yearStr+monStr+dayStr

    # Recently changed this. The hosted data has monthly folders so could add a catch to try that if no file found
    iceTypePathFull = iceTypePathT+'/'+yearStr+'/'+monStr+'/'
    print(iceTypePathFull)
   

    files = glob(iceTypePathFull+'*_nh_*'+dateStr+'1200.nc')
    f = Dataset(files[0], 'r')
    print ('Ice type file path:', files[0])

    lats = f.variables['lat'][::res, ::res]
    lons = f.variables['lon'][::res, ::res]
    ice_typeT = f.variables['ice_type'][0,::res, ::res]
    # set multiyear ice and amiguous to 1, everything else to 0
    ice_typeT = np.where(ice_typeT>2.5, 1, 0)
    ice_typeT = ma.masked_where(lats>88, ice_typeT)

    xpts_type, ypts_type = mapProj(lons, lats)

    return xpts_type, ypts_type, ice_typeT

def getIceType(dF, iceTypePathT, mapProj, res=1, returnRaw=1):
    """
    Get ice type based on the OSI-SAF ice type product, 
    extended back tot he SSMI period as part of the Copernicus Climate Change (C3S) service.
    Thanks to Thomas Lavergne and Singe Aaboe 
    
    Data from 2006 onwards is from OSI-SAF (Lavergene et al.)
    Data prior to 2006 is produced through funding provided by the Copernicus CDR (Aaboe et al.,)
    Product User Guide and Specification (PUGS): Sea Ice Edge and Type version 1, Algorithm Theoretical Basis Document (ATBD): Sea Ice Edge and Type version 2";
    
    Assumes minimal change in ice type across a day so just use the first row for date information and apply to all rows
    
    ice_typeT2 (var): Ice type 

    Args:
        dF (data frame): original data frame
        mapProj (basemap instance): basemap map projection
        res (int): every res number of points sampled along ice type array.
        
    Returns:
        dF (data frame): data frame including ice type column (1 = multiyear ice, 0 = everything else)

    """

    dayStr='%02d' %(dF['day'].iloc[0])
    monthStr='%02d' %(dF['month'].iloc[0])
    yearStr='%02d' %(dF['year'].iloc[0])
    dateStr=yearStr+monthStr+dayStr

    # Recently changed this. THe hosted data has monthly folders so could add a catch to try that if no file found
    iceTypePathFull = iceTypePathT+'/'+yearStr+'/'+monthStr+'/'
    print(iceTypePathFull)
    #if (dF['year'].iloc[0]<2006):
    #    iceTypePathFull = iceTypePathT+'/OSISAF/'+yearStr+'/'
    #    
    #else:
    #    iceTypePathFull = iceTypePathT+'/OSISAF/'+yearStr+'/'+monthStr+'/'
    #print(iceTypePathFull+'*_nh_*'+dateStr+'1200.nc')

    files = glob(iceTypePathFull+'*_nh_*'+dateStr+'1200.nc')
    f = Dataset(files[0], 'r')
    print ('Ice type file path:', files[0])

    lats = f.variables['lat'][::res, ::res]
    lons = f.variables['lon'][::res, ::res]
    ice_typeT = f.variables['ice_type'][0,::res, ::res]
    # set multiyear ice and amiguous to 1, everything else to 0
    ice_typeT = np.where(ice_typeT>2.5, 1, 0)

    xpts_type, ypts_type = mapProj(lons, lats)
    xpts_type=xpts_type.flatten()
    ypts_type=ypts_type.flatten()
    
    ice_typeT2=ice_typeT.flatten()

    #iceTypeGs=[]
    iceTypeGs=ma.masked_all((size(dF['freeboard'].values)))
    for x in range(size(dF['freeboard'].values)):
        # Find nearest ice type
        dist=sqrt((xpts_type-dF['xpts'].iloc[x])**2+(ypts_type-dF['ypts'].iloc[x])**2)
        index_min = np.argmin(dist)
        iceTypeGs[x]=ice_typeT2[index_min]
        # This is waht I normally do but it's much slower!
        # I checked and they gave the same answers
        # iceTypeG2 = griddata((xpts_type, ypts_type), ice_typeT2, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        # print(iceTypeG)
        # iceTypeGs.append(iceTypeG)
    print('good')
    dF['ice_type'] = pd.Series(iceTypeGs, index=dF.index)

    if (returnRaw==1):
        return dF, lats, lons, ice_typeT
    else:
        return dF

def getIceTypeA(dF, iceTypePathT, mapProj, res=1, returnRaw=1):
    """
    Get ice type based on the AMSRE ice type product.
    Available 
    
    Data from 2006 onwards is from OSI-SAF (Lavergene et al.)
    Data prior to 2006 is produced through funding provided by the Copernicus CDR (Aaboe et al.,)
    Product User Guide and Specification (PUGS): Sea Ice Edge and Type version 1, Algorithm Theoretical Basis Document (ATBD): Sea Ice Edge and Type version 2";
    
    Assumes minimal change in ice type across a day so just use the first row for date information and apply to all rows
    
    ice_typeT2 (var): Ice type 

    Args:
        dF (data frame): original data frame
        mapProj (basemap instance): basemap map projection
        res (int): every res number of points sampled along ice type array.
        
    Returns:
        dF (data frame): data frame including ice type column (1 = multiyear ice, 0 = everything else)

    """

    dayStr='%02d' %(dF['day'].iloc[0])
    monthStr='%02d' %(dF['month'].iloc[0])
    yearStr='%02d' %(dF['year'].iloc[0])
    dateStr=yearStr+monthStr+dayStr

    # Recently changed this. THe hosted data has monthly folders so could add a catch to try that if no file found
    iceTypePathFull = iceTypePathT+'/'+yearStr+'/'+monthStr+'/'
    print(iceTypePathFull)
    #if (dF['year'].iloc[0]<2006):
    #    iceTypePathFull = iceTypePathT+'/OSISAF/'+yearStr+'/'
    #    
    #else:
    #    iceTypePathFull = iceTypePathT+'/OSISAF/'+yearStr+'/'+monthStr+'/'
    #print(iceTypePathFull+'*_nh_*'+dateStr+'1200.nc')

    files = glob(iceTypePathFull+'*_nh_*'+dateStr+'1200.nc')
    f = Dataset(files[0], 'r')
    print ('Ice type file path:', files[0])

    lats = f.variables['lat'][::res, ::res]
    lons = f.variables['lon'][::res, ::res]
    ice_typeT = f.variables['ice_type'][0,::res, ::res]
    # set multiyear ice and amiguous to 1, everything else to 0
    ice_typeT = np.where(ice_typeT>2.5, 1, 0)

    xpts_type, ypts_type = mapProj(lons, lats)
    xpts_type=xpts_type.flatten()
    ypts_type=ypts_type.flatten()
    
    ice_typeT2=ice_typeT.flatten()

    #iceTypeGs=[]
    iceTypeGs=ma.masked_all((size(dF['freeboard'].values)))
    for x in range(size(dF['freeboard'].values)):
        # Find nearest ice type
        dist=sqrt((xpts_type-dF['xpts'].iloc[x])**2+(ypts_type-dF['ypts'].iloc[x])**2)
        index_min = np.argmin(dist)
        iceTypeGs[x]=ice_typeT2[index_min]
        # This is waht I normally do but it's much slower!
        # I checked and they gave the same answers
        # iceTypeG2 = griddata((xpts_type, ypts_type), ice_typeT2, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        # print(iceTypeG)
        # iceTypeGs.append(iceTypeG)
    print('good')
    dF['ice_type'] = pd.Series(iceTypeGs, index=dF.index)

    if (returnRaw==1):
        return dF, lats, lons, ice_typeT
    else:
        return dF

#@jit
def getWarrenData(dF, outSnowVar, outDensityVar='None', modFactor=1):
    """
    Assign Warren1999 snow dept/density climatology to dataframe

    Added 

    Args:
        dF (data frame): Pandas dataframe
        outSnowVar (string): name of Warren snow depth variable
        outDensityVar (string): name of Warren snow density variable
        

    Returns:
        dF (data frame): Pandas dataframe updated to include colocated Warren snow depth and density
        
    """
     # 2   Sea of Okhotsk and Japan
    # 3   Bering Sea
    # 4   Hudson Bay
    # 5   Gulf of St. Lawrence
    # 6   Baffin Bay/Davis Strait/Labrador Sea
    # 7   Greenland Sea
    # 8   Barents Seas
    # 9   Kara
    # 10   Laptev
    # 11   E. Siberian
    # 12   Chukchi
    # 13   Beaufort
    # 14   Canadian Archipelago
    # 15   Arctic Ocean

    
    good_regions=[9, 10, 11, 12, 13, 15]

    # Generate empty lists
    snowDepthW99s=ma.masked_all(size(dF['freeboard'].values))
    if (outDensityVar!='None'):
        snowDensityW99s=ma.masked_all(size(dF['freeboard'].values))

    # Loop over all freeboard values (rows)
    for x in range(size(dF['freeboard'].values)):
        #print(x, dF['lon'].iloc[x], dF['lat'].iloc[x], dF['month'].iloc[x]-1)
       
        
        if (dF['region_flag'].iloc[x] in good_regions):

            # Subtract 1 from month as warren index in fucntion starts at 0
            snowDepthDayW99T, snowDensityW99T=WarrenClimatology(dF['lon'].iloc[x], dF['lat'].iloc[x], dF['month'].iloc[x]-1)
            #print(snowDepthDayW99T, snowDensityW99T) 
            # If modFactor less than 1 then this implies we want to modify the climatology based on ice type
            if (modFactor<0.9):
                # If first year ice (should be zero but do less than 0.5 to be safe)
                #print('modify')
                if (dF['ice_type'].iloc[x]<0.5):
                    #print('FYI')
                    #print ('Unmodified snow:', snowDepthDayW99T)
                    snowDepthDayW99T=snowDepthDayW99T*modFactor
                    #print ('Modified snow:', snowDepthDayW99T)

            # Append values to list
            snowDepthW99s[x]=snowDepthDayW99T
            if (outDensityVar!='None'):
                snowDensityW99s[x]=snowDensityW99T
        else:
            # If outside the good regions, apply nan snow depth/density
            snowDepthW99s[x]=np.nan
            if (outDensityVar!='None'):
                snowDensityW99s[x]=np.nan

    # Assign list to dataframe as a series
    dF[outSnowVar] = pd.Series(snowDepthW99s, index=dF.index)
    if (outDensityVar!='None'):
        dF[outDensityVar] = pd.Series(snowDensityW99s, index=dF.index)
    

    return dF

def getWarrenDataCPOM(dF, outSnowVar='W99mod5r', outDensityVar='W99r', modFactor=0.5):
    """
    Assign mean ice type, monthly Warren1999 snow depth/density climatology to dataframe

    Added 

    Args:
        dF (data frame): Pandas dataframe
        outSnowVar (string): name of Warren snow depth variable
        outDensityVar (string): name of Warren snow density variable

    Returns:
        dF (data frame): Pandas dataframe updated to include colocated Warren snow depth and density
        
    """

    # Generate empty lists
    snowDepthW99s=ma.masked_all(size(dF['freeboard'].values))
    if (outDensityVar!='None'):
        snowDensityW99s=ma.masked_all(size(dF['freeboard'].values))

    # Loop over all freeboard values (rows)
    for x in range(size(dF['freeboard'].values)):
        #print(x, dF['lon'].iloc[x], dF['lat'].iloc[x], dF['month'].iloc[x]-1)
       
        snowDepthDayW99T, snowDensityW99T=MonthlyMeanWarrenClimatology(dF['month'].iloc[x]-1)
        
        if (modFactor<0.9):
            # If first year ice (should be zero but do less than 0.5 to be safe)
            #print('modify')
            if (dF['ice_type'].iloc[x]<0.5):
                #print('FYI')
                #print ('Unmodified snow:', snowDepthDayW99T)
                snowDepthDayW99T=snowDepthDayW99T*modFactor
                #print ('Modified snow:', snowDepthDayW99T)


        # Append values to list
        snowDepthW99s[x]=snowDepthDayW99T
        if (outDensityVar!='None'):
            snowDensityW99s[x]=snowDensityW99T
    
    # Assign list to dataframe as a series
    dF[outSnowVar] = pd.Series(snowDepthW99s, index=dF.index)
    if (outDensityVar!='None'):
        dF[outDensityVar] = pd.Series(snowDensityW99s, index=dF.index)
    

    return dF

def MonthlyMeanWarrenClimatology(monthT):
    """
    Get Warren1999 snow depth climatology

    Args:
        lonT (var): longitude
        latT (var): latitude
        monthT (var): month with the index starting at 0
        

    Returns:
        Hs (var): Snow depth (m)
        rho_s (var): Snow density (kg/m^3)
        
    """

    meanCentralArcticSnow = [26.7, 28.9, 31.6, 31.3, 32.6, 27.7, 3.4, 2.5, 8.4, 17.2, 21.4, 23.3]
    meanCentralArcticDensity = [287.5, 305.2, 315.4, 309.8, 321.0, 332.0, 358.8, 238.3, 236.0, 268.2, 280.0, 288.2]

    #meanCentralArcticSnow = [27.1, 29.7, 32.4, 33.7, 34.4, 30.7, 6.6, 3.0, 11.2, 19.4, 23.0, 24.9]
    #meanCentralArcticSWE = [8.0, 9.2, 10.3, 10.6, 11.0, 10.4, 2.4, 0.7, 2.7, 5.2, 6.6, 7.3]

    rho_s = meanCentralArcticDensity[monthT]

    # Density in kg/m^3
    #rho_s = 1000.*(swe/Hs)  
    #print(ma.mean(rho_s))


    Hs = meanCentralArcticSnow[monthT]
    

    # Convert snow depth to meters
    Hs=Hs/100.

    return Hs, rho_s

def WarrenClimatology(lonT, latT, monthT):
    """
    Get Warren1999 snow depth climatology

    Args:
        lonT (var): longitude
        latT (var): latitude
        monthT (var): month with the index starting at 0
        

    Returns:
        Hs (var): Snow depth (m)
        rho_s (var): Snow density (kg/m^3)
        
    """

    H_0 = [28.01, 30.28, 33.89, 36.8, 36.93, 36.59, 11.02, 4.64, 15.81, 22.66, 25.57, 26.67]
    a = [.127, .1056, .5486, .4046, .0214, .7021, .3008, .31, .2119, .3594, .1496, -0.1876]
    b = [-1.1833, -0.5908, -0.1996, -0.4005, -1.1795, -1.4819, -1.2591, -0.635, -1.0292, -1.3483, -1.4643, -1.4229]
    c = [-0.1164, -0.0263, 0.0280, 0.0256, -0.1076, -0.1195, -0.0811, -0.0655, -0.0868, -0.1063, -0.1409, -0.1413]
    d = [-0.0051, -0.0049, 0.0216, 0.0024, -0.0244, -0.0009, -0.0043, 0.0059, -0.0177, 0.0051, -0.0079, -0.0316]
    e = [0.0243, 0.0044, -0.0176, -0.0641, -0.0142, -0.0603, -0.0959, -0.0005, -0.0723, -0.0577, -0.0258, -0.0029]

    # Convert lat and lon into degrees of arc, +x axis along 0 degrees longitude and +y axis along 90E longitude
    x = (90.0 - latT)*cos(lonT * pi/180.0)  
    y = (90.0 - latT)*sin(lonT*pi/180.0) 

    Hs = H_0[monthT] + a[monthT]*x + b[monthT]*y + c[monthT]*x*y + (d[monthT]*x*x) + (e[monthT]*y*y)
    

    # Now get SWE, although this is not returned by the function

    H_0swe = [8.37, 9.43,10.74,11.67,11.8,12.48,4.01,1.08,3.84,6.24,7.54,8.0]
    aswe = [-0.027,0.0058,0.1618,0.0841,-0.0043,0.2084,0.097,0.0712,0.0393,0.1158,0.0567,-0.054]
    bswe = [-0.34,-0.1309,0.0276,-0.1328,-0.4284,-0.5739,-0.493,-0.145,-0.2107,-0.2803,-0.3201,-0.365]
    cswe = [-0.0319,0.0017,0.0213,0.0081,-0.038,-0.0468,-0.0333,-0.0155,-0.0182,-0.0215,-0.0284,-0.0362]
    dswe = [-0.0056,-0.0021,0.0076,-0.0003,-0.0071,-0.0023,-0.0026,0.0014,-0.0053,0.0015,-0.0032,-0.0112]
    eswe = [-0.0005,-0.0072,-0.0125,-0.0301,-0.0063,-0.0253,-0.0343,0,-0.019,-0.0176,-0.0129,-0.0035]


    swe = H_0swe[monthT] + aswe[monthT]*x + bswe[monthT]*y + cswe[monthT]*x*y + dswe[monthT]*x*x + eswe[monthT]*y*y

    # Density in kg/m^3
    rho_s = 1000.*(swe/Hs)  
    #print(ma.mean(rho_s))

    # Could mask out bad regions (i.e. land) here if desired.
    # Hsw[where(region_maskG<9.6)]=np.nan
    # Hsw[where(region_maskG==14)]=np.nan
    # Hsw[where(region_maskG>15.5)]=np.nan

    # Could mask out bad regions (i.e. land) here if desired.
    #rho_s[where(region_maskG<9.6)]=np.nan
    #rho_s[where(region_maskG==14)]=np.nan
    #rho_s[where(region_maskG>15.5)]=np.nan

    # Convert snow depth to meters
    Hs=Hs/100.

    return Hs, rho_s

def getDate(year, month, day):
    """ Get date string from year month and day"""

    return str(year)+'%02d' %month+'%02d' %day


def freeboard_to_thickness(freeboardT, snow_depthT, snow_densityT, ice_densityT):
    """
    Hydrostatic equilibrium equation to calculate sea ice thickness 
    from freeboard and snow depth/density data

    Args:
        freeboardT (var): ice freeboard
        snow_depthT (var): snow depth
        snow_densityT (var): final snow density
        ice_densityT (int): sea ice density

    Returns:
        ice_thicknessT (var): ice thickness dereived using hydrostatic equilibrium

    """

    # Define density values
    rho_w=1024.

    # set snow to freeboard where it's bigger than freeboard.
    snow_depthT[snow_depthT>freeboardT]=freeboardT[snow_depthT>freeboardT]

    ice_thicknessT = (rho_w/(rho_w-ice_densityT))*freeboardT - ((rho_w-snow_densityT)/(rho_w-ice_densityT))*snow_depthT

    return ice_thicknessT

def getIS1FreeboardData(freeboardFileT, mapProj):
    """
    Load ICESat reeboard data into a Pandas dataframe

    Args:
        freeboardFileT (file): file path of ICESat freeboard data
        mapProj (basemap instance): basemap map projection
    Returns:
        dF (var): Dataframe containing freeboard, year, month, day, lon, lat, x, y

    I think the dates are indexed starting from 1 - i.e. month of 1 = January 
        
    """
    

    ds = xr.open_dataset(freeboardFileT)
    dF = ds.to_dataframe()

    dF['year'] = dF.year.astype(int)
    dF['month'] = dF.month.astype(int)
    dF['day'] = dF.day.astype(int)

    # Filter out these negative (NaN) values
    # NB: should we also filter out freeboards equal to zero?
    
    #dF = dF[(dF['freeboard']>0)]
    
    # Set negative values to 0 like IS2
    dF[(dF['freeboard']<0)]=0
    
    # Reset row indexing
    dF=dF.reset_index(drop=True)

    xpts, ypts=mapProj(dF['lon'].values, dF['lat'].values)

    dF['xpts'] = pd.Series(xpts, index=dF.index)
    dF['ypts'] = pd.Series(ypts, index=dF.index)

    #dF2 = dask.dataframe.from_pandas(dF, chunks=1000)
    #dF2.head()

    return dF




def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    Args:
        lon1 (var): initial longitude
        lat1 (var): initial latitude
        lon2 (var): final longitude
        lon2 (var): final longitude

    Returns:
        distance (var): great circle distance in kilometers


    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371 * c
    return km



def snowDistribution(snowDepthT, freeboardT, meanFreeboardT):
    """
    Snow distribution function - simplified version of what Nathan created for testing.
    Written by Alek (summer 2018)

    Args:
        snowDepthT (var): along track snow depth (m)
        freeboardT (var): along track freeboard (m)
        meanFreeboardT (var): mean freebaord (m)

    Returns:
        snowDepthDist (var): along track snow depth distributed to higher resolution
    """

    #print ('mean freeboard section:', meanFreeboardT)

    hs_thick_adj=0.
    fb_cutoff = 0.69*100*snowDepthT +0.22*100*meanFreeboardT + 19.6
    hs_thick = 1.03*snowDepthT*100 + 0.83 + hs_thick_adj
    hs_thin = hs_thick * meanFreeboardT*100.0/fb_cutoff

    if (freeboardT > fb_cutoff):
        snowDepthDist = hs_thick
    else:
        snowDepthDist = hs_thin

    return snowDepthDist/100.

def get_region_mask(ancDataPath, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = ancDataPath+'/region_n.msk'
    
    #1 Non-regional ocean  
    #2 Sea of Okhotsk 
    #3 Bering Sea  
    #4 Hudson Bay 
    #5 Baffin Bay/Davis Strait/Labrador Sea    
    #6 Greenland Sea   Bellingshausen 
    #7 Kara and Barents Seas

	#8 - Arctic Ocean
	#9 - Canadian Archipelago
	#10 - Gulf of St Lawrence
	#11 - Land

	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	if (xypts_return==1):
		mask_latf = open(ancDataPath+'/psn25lats_v3.dat', 'rb')
		mask_lonf = open(ancDataPath+'/psn25lons_v3.dat', 'rb')
		lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask

def get_region_maskCAsnowIS(datapathT, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = datapathT+'/OTHER/region_n.msk'
	
	region_lonlat = [-140, -10, 76, 87]
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapathT+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapathT+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskCA=np.zeros((lons_mask.shape))
	mask = where((lons_mask>region_lonlat[0]) & (lons_mask<region_lonlat[1]) & (lats_mask>region_lonlat[2]) & (lats_mask<region_lonlat[3])& (region_mask==8))
	region_maskCA[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskCA, xpts, ypts
	else:
		return region_maskCA

#@jit
def bindataN(x, y, z, xG, yG, binsize, retbin=True, retloc=False):
    """
    Place unevenly spaced 2D data on a grid by 2D binning (nearest
    neighbor interpolation).

    MAKE SURE YOU ARE USING THE NATIVE PROJECTION OF THE 2D DATA 
    SO THAT IT IS ORDERD CORRECTLY TO GENERATE THE MESHGRID BELOW.
    
    Parameters
    ----------
    x : ndarray (1D)
        The idependent data x-axis of the grid.
    y : ndarray (1D)
        The idependent data y-axis of the grid.
    z : ndarray (1D)
        The dependent data in the form z = f(x,y).
    binsize : scalar, optional
        The full width and height of each bin on the grid.  If each
        bin is a cube, then this is the x and y dimension.  This is
        the step in both directions, x and y. Defaults to 0.01.
    retbin : boolean, optional
        Function returns `bins` variable (see below for description)
        if set to True.  Defaults to True.
    retloc : boolean, optional
        Function returns `wherebins` variable (see below for description)
        if set to True.  Defaults to True.
   
    Returns
    -------
    grid : ndarray (2D)
        The evenly gridded data.  The value of each cell is the median
        value of the contents of the bin.
    bins : ndarray (2D)
        A grid the same shape as `grid`, except the value of each cell
        is the number of points in that bin.  Returns only if
        `retbin` is set to True.
    wherebin : list (2D)
        A 2D list the same shape as `grid` and `bins` where each cell
        contains the indicies of `z` which contain the values stored
        in the particular bin.

    Revisions
    ---------
    2010-07-11  ccampo  Initial version
    """
    # get extrema values.
    xmin, xmax = xG.min(), xG.max()
    ymin, ymax = yG.min(), yG.max()

    # make coordinate arrays.
    xi      = xG[0]
    yi      = yG[:, 0] #np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)
    print(xi)
    print(yi)
    # make the grid.
    grid           = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)

    # create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()

    # fill in the grid.
    for row in range(nrow):
        for col in range(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.

            # find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
            ind  = np.where(ibin == True)[0]

            # fill the bin.
            bin = z[ibin]
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size != 0:
                binval         = np.mean(bin)
                grid[row, col] = binval
            else:
                grid[row, col] = np.nan   # fill empty bins with nans.

    # return the grid
    if retbin:
        if retloc:
            return grid, bins, wherebin
        else:
            return grid, bins
    else:
        if retloc:
            return grid, wherebin
        else:
            return grid

#@jit
def bindataSegWeighted(x, y, z, seg, xG, yG, binsize=0.01, retbin=True, retloc=False):
    """
    Place unevenly spaced 2D data on a grid by 2D binning (nearest
    neighbor interpolation) and weight using the IS2 segment lengths.
    
    Parameters
    ----------
    x : ndarray (1D)
        The idependent data x-axis of the grid.
    y : ndarray (1D)
        The idependent data y-axis of the grid.
    z : ndarray (1D)
        The dependent data in the form z = f(x,y).
    seg : ndarray (1D)
        The segment length of the data points in the form z = seg(x,y).
    binsize : scalar, optional
        The full width and height of each bin on the grid.  If each
        bin is a cube, then this is the x and y dimension.  This is
        the step in both directions, x and y. Defaults to 0.01.
    retbin : boolean, optional
        Function returns `bins` variable (see below for description)
        if set to True.  Defaults to True.
    retloc : boolean, optional
        Function returns `wherebins` variable (see below for description)
        if set to True.  Defaults to True.
   
    Returns
    -------
    grid : ndarray (2D)
        The evenly gridded data.  The value of each cell is the median
        value of the contents of the bin.
    bins : ndarray (2D)
        A grid the same shape as `grid`, except the value of each cell
        is the number of points in that bin.  Returns only if
        `retbin` is set to True.
    wherebin : list (2D)
        A 2D list the same shape as `grid` and `bins` where each cell
        contains the indicies of `z` which contain the values stored
        in the particular bin.

    Revisions
    ---------
    2010-07-11  ccampo  Initial version
    """
    # get extrema values.
    xmin, xmax = xG.min(), xG.max()
    ymin, ymax = yG.min(), yG.max()

    # make coordinate arrays.
    xi      = xG[0]
    yi      = yG[:, 0] #np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)

    # make the grid.
    grid           = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)

    # create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()

    # fill in the grid.
    for row in prange(nrow):
        for col in prange(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.

            # find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
            ind  = np.where(ibin == True)[0]

            # fill the bin.
            bin = z[ibin]
            segbin = seg[ibin]
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size != 0:
                binvalseg         = np.sum(bin*segbin)/np.sum(segbin)
                grid[row, col] = binvalseg
            else:
                grid[row, col] = np.nan   # fill empty bins with nans.

    # return the grid
    if retbin:
        if retloc:
            return grid, bins, wherebin
        else:
            return grid, bins
    else:
        if retloc:
            return grid, wherebin
        else:
            return grid

def monLabels(month):
    labels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 
    return labels[month]

def defGrid(m, dxRes=50000):
    # took out the +1 to shorten array slightly
    nx = int((m.xmax-m.xmin)/dxRes); ny = int((m.ymax-m.ymin)/dxRes)
    gridStr=str(int(dxRes/1000))+'km'
    lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)

    return lonsG, latsG, xptsG, yptsG, nx, ny


def read_icebridgeALL(dataPath, year, mask_hi=1, mask_nonarctic=1):
    # Now using Pandas for this (see OIB notebook) which is a lot easier..
    # Good to comapre with this more laborious approach
    # Will add more funcs in due course
    lats_total=[] 
    lons_total=[]
    thickness_total=[]
    snow_thickness_total=[]
    
    files = glob(rawdatapath+'/ICEBRIDGE/ICEBRIDGE_HI/IDCSI4/'+str(year)+'.*/*.txt')
    
    for x in xrange(size(files)):
        data = genfromtxt(files[x], delimiter=',', skip_header=1, dtype=str)
        # data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
        lats = data[:, 0].astype(float)
        lons = data[:, 1].astype(float)
        thickness = data[:, 2].astype(float)
        snow_thickness = data[:, 7].astype(float)
        lats_total.extend(lats)
        lons_total.extend(lons)
        thickness_total.extend(thickness)
        snow_thickness_total.extend(snow_thickness)

    thickness_total=array(thickness_total)
    snow_thickness_total=array(snow_thickness_total)
    lats_total=array(lats_total)
    lons_total=array(lons_total)

    if (mask_hi==1):
        good_data=where((thickness_total>=0.)&(thickness_total<=20.))
        thickness_total = thickness_total[good_data]
        snow_thickness_total=snow_thickness_total[good_data]
        lats_total = lats_total[good_data]
        lons_total = lons_total[good_data]
    if (mask_nonarctic==1):
        region_mask, xptsM, yptsM = get_region_mask(rawdatapath, mplot)
        xptsIB, yptsIB = mplot(lons_total, lats_total)
        region_maskR = griddatascipy((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsIB, yptsIB), method='nearest')
        good_data = where((region_maskR==8))
        lats_total = lats_total[good_data]
        lons_total=lons_total[good_data]
        thickness_total=thickness_total[good_data]
        snow_thickness_total=snow_thickness_total[good_data]

    xpts,ypts = mplot(lons_total, lats_total)

    return xpts,ypts, lats_total, lons_total, thickness_total, snow_thickness_total


def mapll(alat, along):
    # AUthor: Nathan
    # Converts lat and lon into x,y


    # Conversion constant from degrees to radians = 57.29577951.      
    cdr=57.29577951 
    # Radius of Earth = 6378.273 km (Hughes Ellipsoid).  
    re=6378.273 
    # Eccentricity of the Earth   e**2=0.006693883 (Hughes Ellipsoid). 
    e2=0.006693883 
    e=math.sqrt(e2) 

    pi=3.14159265359 
    # Standard Parallel - Latitude with No distortion =  +_ 70 degrees 
    slat=70.0 

    #**** Degree to Radiants :
    d2r=3.14159265359/180.

    if alat > 0.0: # Northern Hemisphere
        xlam=-45.
        sn=1.0
    else:
        xlam=0.
        sn=-1.0


    #     Compute X and Y in grid coordinates. 
    alat=sn*alat 
    along = sn*along


    rlat = alat


    t0=math.tan((pi/4.0)-(rlat/(2.0*cdr)))/((1.0-e*math.sin(rlat*d2r)) /  (1.0+e*math.sin(rlat*d2r)))**(e/2.0) 
    rlat = slat
    t1=math.tan((pi/4.0)-(rlat/(2.0*cdr)))/((1.0-e*math.sin(rlat*d2r))/  (1.0+e*math.sin(rlat*d2r)))**(e/2.0) 
    cm=math.cos(slat*d2r)/math.sqrt(1.0-e2*(math.sin(slat*d2r)**2)) 
    rho=re*cm*t0/t1 
    x=rho*sn*math.sin((along-xlam)*d2r) 
    y=-1.0*rho*sn*math.cos((along-xlam)*d2r) 




    if alat >= 89.995:
        x=0.0
        y=0.0

