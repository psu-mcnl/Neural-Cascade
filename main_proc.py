import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy import linalg
from scipy import stats
from scipy.fftpack import fft, ifft
import scipy.io as sio
import os
import shutil
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import pickle
import copy
import random
from pathlib import Path

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession, 
    removed_unused_stimulus_presentation_columns
)
from allensdk.brain_observatory.ecephys.visualization import plot_mean_waveforms, plot_spike_counts, raster_plot
from allensdk.brain_observatory.visualization import plot_running_speed

from utils import *
from settings import settings
from hyperparams import *

def extract_seg(sig,centers,hwinsz):
    """ extract segments from signals """
    
    if sig.ndim == 2:
        seg = np.zeros((hwinsz*2+1,sig.shape[1],centers.size))
        for cnt in range(centers.size):
            seg[:,:,cnt] = sig[centers[cnt]-hwinsz:centers[cnt]+hwinsz+1,:]
        return seg
    elif sig.ndim == 1:     
        seg = np.zeros((hwinsz*2+1,centers.size))
        for cnt in range(centers.size):
            seg[:,cnt] = sig[centers[cnt]-hwinsz:centers[cnt]+hwinsz+1]
        return seg
    else:
        sys.exit('ndim of sig musts be either 1 or 2')

def delayprofile_centroid(sigorg, binEdgeIdx):
    """ function to compute the principal delay profile using centroid """
    
    sig = sigorg.copy()
    
    sig_min = sig.min(axis=0)
    idx = sig_min < 0
    sig[:,idx] = sig[:,idx] - sig_min[idx]
        
    dprfs = []
    binpks = []
    for sidx, eidx in zip(binEdgeIdx[:-1], binEdgeIdx[1:]):
        sigseg = sig[sidx:eidx,:]
        t = np.arange(eidx - sidx)
        dprfseg = sigseg.T @ t / sigseg.sum(axis=0)
        peak = np.nanmean(dprfseg) + sidx
        
        dprfs.append(dprfseg)
        binpks.append(peak)
    
    dprfs = np.array(dprfs)
    binpks = np.round(binpks).astype('int')
    
    idx0, idx1 = np.isnan(dprfs).nonzero()
    # dprfs[idx0, idx1] = np.nanmean(dprfs,axis=0)[idx1]
    dprfs[idx0, idx1] = np.nanmean(dprfs,axis=1)[idx0]
    
    return dprfs, binpks

def ripple_detection(sig, sr):
    """ return ripple event mask and ripple peak index 
        input:
            sig: lfp data in the form of N time points and M channels
            sr: sampling rate of the lfp signal
    """
    lowB = 80
    highB = 250
    SDthr = 5
    b1, a1 = signal.butter(2, Wn = [lowB/sr*2, highB/sr*2], btype = "bandpass")
    b2, a2 = signal.butter(2, Wn = (lowB+highB)/2/3.1415926/sr*2, btype = "lowpass")
    ripplemsk = []
    ripplectrind = []
    sig = signal.filtfilt(b1, a1, sig, axis=0,method="gust") 
    for cnt in range(sig.shape[1]):
        org = sig[:,cnt]
        clipped = np.clip(org,org.mean()-org.std()*SDthr,org.mean()+org.std()*SDthr)
        clipped_pw = signal.filtfilt(b2, a2, abs(clipped),axis=0,method='gust')
        org_pw = signal.filtfilt(b2, a2, abs(org),axis=0,method='gust')

        sd5msk = (org_pw>(clipped_pw.mean()+clipped_pw.std()*SDthr)).astype(int)
        sd5msk[[0,-1]]=0
        # remove small ripple mask (<15ms)
        sd5msk_diff = np.diff(sd5msk)
        sd5startind = np.asarray(np.where(sd5msk_diff==1)).squeeze()+1
        sd5endind = np.asarray(np.where(sd5msk_diff==-1)).squeeze()+1
        sd5blksz = sd5endind-sd5startind
        for cnt2 in range(sd5blksz.shape[0]):
            if sd5blksz[cnt2] < 15/1000*sr:
                sd5msk[sd5startind[cnt2]:sd5endind[cnt2]]=0
        # merge ripple mask with a small gap (<15ms)
        sd5msk_diff = np.diff(sd5msk)
        sd5startind = np.asarray(np.where(sd5msk_diff==1)).squeeze()+1
        sd5endind = np.asarray(np.where(sd5msk_diff==-1)).squeeze()+1
        if sd5startind.size>1:
            sd5gapsz = sd5startind[1:]-sd5endind[:-1]
            for cnt3 in range(sd5gapsz.size):
                if sd5gapsz[cnt3] < 15/1000*sr:
                    sd5msk[sd5endind[cnt3]:sd5startind[cnt3+1]]=1
        else:
            sd5msk = np.zeros(sd5msk.shape)
            
        ripplepwpksind, _ = signal.find_peaks(sd5msk*org_pw)
        
        sd2msk = (org_pw>(clipped_pw.mean()+clipped_pw.std()*2)).astype(int)
        sd2msk[[0,-1]]=0
        sd2msk_diff = np.diff(sd2msk)
        sd2startind = np.asarray(np.where(sd2msk_diff==1)).squeeze()+1
        sd2endind = np.asarray(np.where(sd2msk_diff==-1)).squeeze()+1
        for cnt4 in range(sd2startind.size):
            if np.sum((ripplepwpksind-sd2startind[cnt4])*(ripplepwpksind-sd2endind[cnt4])<0)==0:
                sd2msk[sd2startind[cnt4]:sd2endind[cnt4]]=0
        print('processed channel [{:d}|{:d}]'.format(cnt+1, sig.shape[1]))
        ripplemsk.append(sd2msk*org_pw)
        
        if ripplepwpksind.size==0:
            ripplectrind.append(np.nan)
        else:
            hfnpks, _ = signal.find_peaks(-org)
            ripplectrind.append(hfnpks[np.abs(np.subtract.outer(hfnpks, ripplepwpksind)).argmin(0)])
            
    return ripplemsk, ripplectrind    

def process_rsData(session_data, use_all_neurons=True):
    rsData = StimulusSpecificData('spontaneous')

    """ resting-state info """
    stimTable = session_data.get_stimulus_table('spontaneous')

    rsData.stimTable = stimTable
    rsData.stimID = stimTable.index[-2]
    rsData.tstart = stimTable['start_time'].values[-2]
    rsData.tstop = stimTable['stop_time'].values[-2]
    rsData.tdur = stimTable['duration'].values[-2]

    """ neuron selection """
    unitid = session_data.units.index.values
    unitloc = session_data.units['ecephys_structure_acronym']
    if use_all_neurons:
        glsnind = np.ones(unitid.shape[0]).astype(np.bool)
    else:
        THInd = [unitloc == area for area in TH_AREA_NAME]
        THInd = np.array(THInd).any(axis=0)
        HIPInd = [unitloc == area for area in HIP_AREA_NAME]
        HIPInd = np.array(HIPInd).any(axis=0)
        VISInd = unitloc.str.match(VIS_AREA_PREFIX).values

        nPerRegion = min([sum(THInd), sum(HIPInd), sum(VISInd)])
        THIdx = np.random.choice(THInd.nonzero()[0], nPerRegion, replace=False)
        HIPIdx = np.random.choice(HIPInd.nonzero()[0], nPerRegion, replace=False)
        VISIdx = np.random.choice(VISInd.nonzero()[0], nPerRegion, replace=False)

        # neuron selected and will contribute global signal
        glsnind = np.zeros(unitid.shape[0]).astype(np.bool)
        glsnind[THIdx] = True
        glsnind[HIPIdx] = True
        glsnind[VISIdx] = True

    print('number of neurons = {:d}'.format(sum(glsnind)))
    rsData.glsnind = glsnind
    unitid = unitid[rsData.glsnind]
    unitloc = unitloc[rsData.glsnind]

    """ process neural spiking activity data """
    tstep = 0.2

    rsData.tstep = tstep
    rsData.tbins = np.arange(0,int(rsData.tdur)+1,tstep)
    rsData.spikect = session_data.presentationwise_spike_counts(rsData.tbins, rsData.stimID, unitid).squeeze()
    rsData.spikect2 = rsData.spikect / rsData.spikect.mean(axis=0)
    rsData.spikect2f = signal.filtfilt(FILTERS.flt1.b, FILTERS.flt1.a, rsData.spikect2, axis=0, method="gust")
    rsData.spikect3 = (rsData.spikect - rsData.spikect.mean(axis=0)) / rsData.spikect.std(axis=0)

    """ global signal """
    rsData.gls = rsData.spikect2.mean(axis=1)
    rsData.gls2 = rsData.spikect2f.mean(axis=1)
    rsData.glsf = signal.filtfilt(FILTERS.flt1.b, FILTERS.flt1.a, rsData.gls, method="gust") 
    rsData.gls2f = signal.filtfilt(FILTERS.flt1.b, FILTERS.flt1.a, rsData.gls2, method="gust") 

    npks, _ = signal.find_peaks(-rsData.glsf)
    npks2, _ = signal.find_peaks(-rsData.gls2f)
    npks = npks[1:-1]
    npks2 = npks2[1:-1]
    rsData.npks = npks
    rsData.npks2 = npks2

    """ process running data """
    running = session_data.running_speed

    stim_running = running[(running['start_time']>rsData.tstart) & (running['end_time']<rsData.tstop)]
    running_speed = stim_running.velocity.values
    running_time = stim_running.iloc[:,0:2].mean(axis=1) - rsData.tstart
    # 793224716 gets a outlier
    running_speed[running_speed < -90 * running_speed.std()] = 0

    # stationary mask during resting state
    running_speedf = signal.filtfilt(FILTERS.flt2.b, FILTERS.flt2.a, running_speed, axis=0, method="gust") 
    thres = -np.percentile(running_speedf[running_speedf<0], .05)
    qtmskorg = (running_speedf < thres).astype(int)
    qtmskorg = np.concatenate(([0], qtmskorg, [0]))
    qtmskorg_diff = np.diff(qtmskorg)
    qtstartind = np.asarray(np.where(qtmskorg_diff==1)).squeeze() + 1
    qtendind = np.asarray(np.where(qtmskorg_diff==-1)).squeeze() + 1
    qtblksz = qtendind - qtstartind
    for cnt in range(qtblksz.shape[0]):
        if qtblksz[cnt] < 2000:
            qtmskorg[qtstartind[cnt]:qtendind[cnt]]=0
    qtmskorg = qtmskorg[1:-1]

    qtmsk,_,_ = stats.binned_statistic(running_time.values, qtmskorg, statistic='mean', bins=rsData.tbins)
    qtmsk[-1] = qtmsk[-2]
    qtmsk = np.round(qtmsk)

    rsData.running_speed = running_speed
    rsData.running_time = running_time.values
    rsData.running_speedf = running_speedf
    rsData.qtmskorg = qtmskorg
    rsData.qtmsk = qtmsk

    """ process pupil data """
    pupil = session_data.get_pupil_data()
    stim_pupil = pupil[(pupil.index>rsData.tstart) & (pupil.index<rsData.tstop)]
    stim_pupil.index = stim_pupil.index - rsData.tstart

    pupil_size = stim_pupil.loc[:, ['pupil_height', 'pupil_width']].mean(axis=1).values
    pupil_time = stim_pupil.index.values

    rsData.pupil_data = stim_pupil
    rsData.pupil_size = pupil_size
    rsData.pupil_time = pupil_time

    """ delay profiles & principal delay profile  """
    dprfs, glspksind = delayprofile_centroid(rsData.spikect.values, rsData.npks)
    dprfs2, glspksind2 = delayprofile_centroid(rsData.spikect2f, rsData.npks)
    dprfsn = dprfs.T - dprfs.mean(axis=1).T
    dprfsn = dprfsn / dprfsn.std(axis=0)
    dprfs2n = dprfs2.T - dprfs2.mean(axis=1).T
    dprfs2n = dprfs2n / dprfs2n.std(axis=0)

    rsData.dprfs = dprfs
    rsData.glspksind = glspksind
    rsData.dprfs2 = dprfs2
    rsData.glspksind2 = glspksind2
    rsData.dprfsn = dprfsn
    rsData.dprfs2n = dprfs2n
    rsData.glspkstime = rsData.spikect.time_relative_to_stimulus_onset[rsData.glspksind]
    rsData.glspkstime2 = rsData.spikect.time_relative_to_stimulus_onset[rsData.glspksind2]

    U, s, Vh = linalg.svd(rsData.dprfsn)
    U2, s2, Vh2 = linalg.svd(rsData.dprfs2n)
    if Vh[0,:].mean()<0:
        U = -U
        Vh = -Vh
    if Vh2[0,:].mean()<0:
        U2 = -U2
        Vh2 = -Vh2

    rsData.U = U
    rsData.Vh = Vh
    rsData.U2 = U2
    rsData.Vh2 = Vh2

    rsData.dlycor = np.corrcoef(rsData.dprfsn.T, rsData.U[:,0].T)[:-1,-1]
    rsData.dlycor2 = np.corrcoef(rsData.dprfs2n.T, rsData.U2[:,0].T)[:-1,-1]
    rsData.sind =np.argsort(rsData.U[:,0])
    rsData.sind2 =np.argsort(rsData.U2[:,0])

    """ construct unit dataframe and channel dataframe with attribute (delay profile, position) """
    unitdf = units.loc[unitid, ['ecephys_channel_id']]
    unitdf['delays'] = U[:,0]
    unitdf['apc'] = units.loc[unitid, 'anterior_posterior_ccf_coordinate']
    unitdf['dvc'] = units.loc[unitid, 'dorsal_ventral_ccf_coordinate']
    unitdf['lrc'] = units.loc[unitid, 'left_right_ccf_coordinate']
    unitdf['vc'] = units.loc[unitid, 'probe_vertical_position']
    unitdf['loc'] = units.loc[unitid, 'ecephys_structure_acronym']

    chndf = unitdf[['ecephys_channel_id', 'delays']].groupby('ecephys_channel_id').mean()
    chnid = chndf.index.values
    chndf['apc'] = channels.loc[chnid, 'anterior_posterior_ccf_coordinate']
    chndf['dvc'] = channels.loc[chnid, 'dorsal_ventral_ccf_coordinate']
    chndf['lrc'] = channels.loc[chnid, 'left_right_ccf_coordinate']
    chndf['vc'] = channels.loc[chnid, 'probe_vertical_position']
    chndf['loc'] = channels.loc[chnid, 'ecephys_structure_acronym']

    rsData.unitdf = unitdf
    rsData.chndf = chndf

    return rsData

def process_lfpData(session_data, rsData):

    lfpData = DataContainer()

    probesid = session_data.probes.index.values
    sr = session_data.probes.lfp_sampling_rate.values[0]    

    lfp_time = None
    t_pad = 0

    rs_lfpCA1 = list()
    rs_lfpVIS = list()
    for pid in probesid:
        print('processing lfp data for probe-{:d}'.format(pid))
        rs_lfp = session_data.get_lfp(pid).sel(time=slice(rsData.tstart-t_pad, rsData.tstop+t_pad))
        rs_lfp = rs_lfp.assign_coords(time = rs_lfp.time - rsData.tstart)
        
        chnid = rs_lfp.channel.values
        chnloc = channels.loc[chnid, 'ecephys_structure_acronym']

        rs_lfp = rs_lfp.drop_sel(channel=chnid[np.where(chnloc.str.match('NaN'))])
        chnid = rs_lfp.channel.values
        chnloc = channels.loc[chnid, 'ecephys_structure_acronym']
        
        if lfp_time is not None:
            intp = interp1d(rs_lfp.time.values.T, rs_lfp.values.T)
            sig_intp = intp(lfp_time.values).T
            rs_lfp = xr.DataArray(sig_intp, coords=[lfp_time, rs_lfp.channel])
        else:
            lfp_time = rs_lfp.time
            t_pad = 1
        
        rs_lfpCA1.append(rs_lfp.sel(channel=chnid[np.where(chnloc.str.match('CA1'))]))
        rs_lfpVIS.append(rs_lfp.sel(channel=chnid[np.where(chnloc.str.match('VIS'))]))

    rs_lfpCA1 = xr.concat(rs_lfpCA1, dim="channel")
    rs_lfpVIS = xr.concat(rs_lfpVIS, dim="channel")

    rs_lfpCA1 = rs_lfpCA1.fillna(0)
    rs_lfpVIS = rs_lfpVIS.fillna(0)

    lfpData.probesid = probesid
    lfpData.sr = sr

    """ ripple detection using lfp data from CA1 region """
    ripplemsk, ripplectrind = ripple_detection(rs_lfpCA1.values.copy(), sr)

    ripplemskind = []
    for ripple in ripplemsk:
        ripplemskind.append(np.array(np.where(ripple!=0)).squeeze())
    ripplemsk = np.array(ripplemsk).T

    CA1chn = rs_lfpCA1.channel.values
    CA1time = lfp_time.values

    ripplemskall = (ripplemsk > 0).mean(axis=1)
    chn = ripplemsk[ripplemskall > 0.4,:].mean(axis = 0).argmax()

    rprawseg = extract_seg(rs_lfpCA1.values[:,chn], ripplectrind[chn][1:-1], 600)
    rprawseg_t = np.arange(-600, 601) / sr

    lfpData.ripplemsk = ripplemsk
    lfpData.ripplectrind = ripplectrind
    lfpData.ripplemskind = ripplemskind
    lfpData.CA1chn = CA1chn
    lfpData.CA1time = CA1time
    lfpData.chn = chn
    lfpData.ripplemskall = ripplemskall
    lfpData.rprawseg = rprawseg
    lfpData.rprawseg_t = rprawseg_t

    """ get delta power through visual channels """
    rs_lfpVISflt = signal.filtfilt(FILTERS.flt3.b, FILTERS.flt3.a, rs_lfpVIS, axis=0, method="gust") 

    rs_lfpVISdp = signal.filtfilt(FILTERS.flt4.b, FILTERS.flt4.a, abs(rs_lfpVISflt), axis=0, method="gust") 
    rs_lfpVISdp = xr.DataArray(rs_lfpVISdp, coords=[rs_lfpVIS.time, rs_lfpVIS.channel])
    rs_lfpVISdp = rs_lfpVISdp.coarsen(time=10, boundary='trim').mean()

    lfpData.rs_lfpVISdp = rs_lfpVISdp

    return lfpData

def process_stimData(session_data):
    unitid = (session_data.units.index.values)
    unitloc = session_data.units['ecephys_structure_acronym']

    """ stimulus data: flashes """
    stimTable = session_data.get_stimulus_table('flashes')
    flData = StimulusSpecificData('flashes')
    flData.stimTable = stimTable

    tstep = 0.05
    flData.tstep = tstep
    flData.tbins = np.arange(-0.5, 1.5, tstep)
    flData.spikect = session_data.presentationwise_spike_counts(flData.tbins, stimTable.index, unitid)

    """ stimulus data: drifting_gratings_75_repeats """
    stimTable = session_data.get_stimulus_table('drifting_gratings_75_repeats')
    dgrData = StimulusSpecificData('drifting_gratings_75_repeats')
    dgrData.stimTable = stimTable

    tstep = 0.05
    dgrData.tstep = tstep
    dgrData.tbins = np.arange(-1, 2, tstep)
    dgrData.spikect = session_data.presentationwise_spike_counts(dgrData.tbins, stimTable.index, unitid)

    """ stimulus data: drifting_gratings_contrast """
    stimTable = session_data.get_stimulus_table('drifting_gratings_contrast')
    dgcData = StimulusSpecificData('drifting_gratings_contrast')
    dgcData.stimTable = stimTable

    tstep = 0.05
    dgcData.tstep = tstep
    dgcData.tbins = np.arange(-.5, .5, tstep)
    dgcData.spikect = session_data.presentationwise_spike_counts(dgcData.tbins, stimTable.index, unitid)

    return flData, dgrData, dgcData

def process_session(sid, vars2proc):
    """ session data """
    session_data = cache.get_session_data(sid)
    unitid = (session_data.units.index.values)
    unitloc = session_data.units['ecephys_structure_acronym']

    sessData = DataContainer()
    sessData.unitid = unitid
    sessData.unitloc = unitloc

    wvforms = list()
    unit2chn = session_data.units.loc[:, 'peak_channel_id']
    for i, uid in enumerate(unitid):
        wv = session_data.mean_waveforms[uid]
        wv = wv.sel(channel_id=unit2chn[uid])
        wvforms.append(wv)
    wvforms = np.array(wvforms)
    t_wvform = session_data.mean_waveforms[unitid[0]].time.values

    sessData.wvforms = wvforms
    sessData.t_wvform = t_wvform

    if vars2proc.rsData1 | vars2proc.lfpData:
        print('processing resting-state data (all neurons)')
        rsData1 = process_rsData(session_data, use_all_neurons=True)

    if vars2proc.rsData2:
        if sid not in settings.observatory.sess_few_TH_units:
            print('processing resting-state data (evenly distributed)')
            rsData2 = process_rsData(session_data, use_all_neurons=False)

    if vars2proc.flData | vars2proc.dgrData | vars2proc.dgcData:
        print('processing stimulus data')
        flData, dgrData, dgcData = process_stimData(session_data)

    if vars2proc.lfpData:
        print('processing lfp data')
        lfpData = process_lfpData(session_data, rsData1)

    """ save data """
    sessDir = settings.projectData.dir.sessions / 'session_{:d}'.format(sid)
    if not os.path.exists(sessDir):
        print('create directory: {:s}'.format(str(sessDir)))
        os.makedirs(sessDir, exist_ok=True)

    var2save = [key for key in vars2proc.keys() if vars2proc[key]]
    if sid in settings.observatory.sess_few_TH_units:
        var2save.remove('rsData2')

    for varName in var2save:
        fName = getattr(settings.projectData.files.sessions, varName)
        fPath = sessDir / fName

        with open(fPath, 'wb') as handle:
            print('save {:s} to {:s} ...'.format(varName, str(fPath)))
            pickle.dump(eval(varName), handle, protocol=pickle.HIGHEST_PROTOCOL)
            print('... done')

cache = EcephysProjectCache.from_warehouse(manifest=settings.observatory.manifest_path)
sessions = cache.get_session_table()
channels = cache.get_channels()
units = cache.get_units()

session_fc_ids = sessions[sessions['session_type']=='functional_connectivity'].index.values
session_used = np.setdiff1d(session_fc_ids,settings.observatory.sess_removed)

def main():
    vars2proc = AttrDict(
        rsData1 = True,
        flData = False,
        dgrData = False,
        dgcData = False,
        lfpData = False,
        sessData = True,
        rsData2 = True
    )

    NSess = len(session_used)
    for i, sid in enumerate(session_used):
        print('process session [{:d}|{:d}]'.format(i+1, NSess))
        process_session(sid, vars2proc)
    
    # process_session(session_used[1], vars2proc)

if __name__ == "__main__":
    main()
