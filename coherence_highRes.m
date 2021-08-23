clc; clear;

%% general
homeDir = '/home/yifan/';
projectDataDir = [homeDir, 'Data/Project-Data/Neuroscience/Neural-Cascade/'];
sessDir = [projectDataDir, 'sessions/'];

filename = 'rs_spike_data.mat';

session_used = [766640955, 771160300, 771990200, 774875821, 778998620, ...
                779839471, 781842082, 793224716, 794812542, 816200189, ...
                819186360, 821695405, 829720705, 847657808];
            
%% session loop           
for si =1:length(session_used)
    sid = session_used(si);
    filepath = [sessDir, sprintf('session_%d/', sid), filename];
    vars = load(filepath);
    
    nNeurons = length(vars.spiket);
    
    %% coherence between binned point processes (stationary, binned, low-res)
    USE_FILTERED_SIGNAL = false;
    fprintf('[%d] compute coherence (stationary, binned, low-res)\n', sid)
    
    clear params cohpb1
    params.Fs = 50; % sampling frequency
    params.fpass = [0 25]; % frequency range of interest
    params.tapers = [5 9]; % tapers
    params.err = [2 0.05];
    params.pad = 0;
    win = 75;
    
    m = logical(vars.rsData1_qtmsk);
    if sum(m)/params.Fs < win
        win = sum(m)/params.Fs;
    end
    cohall = [];
    phiall = [];
    for i = 1:nNeurons
        fprintf('processing [%d|%d]\n', i, nNeurons);
        if USE_FILTERED_SIGNAL
            spikect = vars.rsData1_spikect2f(m,:);
        else
            spikect = vars.rsData1_spikect2(m,:);
        end
        if size(spikect,1) > 10000
            spikect = spikect(2500:end-2500,:);
        end
        sig = spikect(:,i);
        gls = [spikect(:,1:i-1) spikect(:,i+1:end)];
        gls = nanmean(gls, 2);
        
        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr] = coherencysegpb(sig, gls, win, params, 1);
        
        cohall = [cohall C];
        phiall = [phiall phi];    
    end
    
    cohpb1.cohall = cohall;
    cohpb1.phiall = phiall;
    cohpb1.f = f;
    
    %% coherence between binned point processes (stationary, binned, low-res, phase-scrambled)
    fprintf('[%d] compute coherence (stationary, binned, low-res, phase-scrambled)\n', sid)
    
    clear params cohpb2 m
    params.Fs = 50; % sampling frequency
    params.fpass = [0 25]; % frequency range of interest
    params.tapers = [5 9]; % tapers
    params.err = [2 0.05];
    params.pad = 0;
    win = 75;
    
    if size(vars.rsData1_spikect2_qtps,1)/params.Fs < win
        win = size(vars.rsData1_spikect2_qtps,1)/params.Fs;
    end
    cohall = [];
    phiall = [];
    for i = 1:nNeurons
        fprintf('processing [%d|%d]\n', i, nNeurons);        
        spikect = vars.rsData1_spikect2_qtps;
        if size(spikect,1) > 10000
            spikect = spikect(2500:end-2500,:);
        end
        
        sig = spikect(:,i);
        gls = [spikect(:,1:i-1) spikect(:,i+1:end)];
        gls = nanmean(gls, 2);
    
        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr] = coherencysegpb(sig, gls, win, params, 1);
        cohall = [cohall C];
        phiall = [phiall phi];    
    end
    
    cohpb2.cohall = cohall;
    cohpb2.phiall = phiall;
    cohpb2.f = f;

    %% coherence between binned point processes (running, binned, low-res)
    USE_FILTERED_SIGNAL = false;
    fprintf('[%d] compute coherence (running, binned, low-res)\n', sid)

    clear params cohpb3
    params.Fs = 50; % sampling frequency
    params.fpass = [0 25]; % frequency range of interest
    params.tapers = [5 9]; % tapers
    params.err = [2 0.05];
    params.pad = 0;
    win = 75;
    
    m = ~logical(vars.rsData1_qtmsk);
    if sum(m)/params.Fs < win
        win = sum(m)/params.Fs;
    end
    cohall = [];
    phiall = [];
    for i = 1:nNeurons
        fprintf('processing [%d|%d]\n', i, nNeurons);
        if USE_FILTERED_SIGNAL
            spikect = vars.rsData1_spikect2f(m,:);
        else
            spikect = vars.rsData1_spikect2(m,:);
        end
        if size(spikect,1) > 10000
            spikect = spikect(2500:end-2500,:);
        end
        sig = spikect(:,i);
        gls = [spikect(:,1:i-1) spikect(:,i+1:end)];
        gls = nanmean(gls, 2);
        
        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr] = coherencysegpb(sig, gls, win, params, 1);
        
        cohall = [cohall C];
        phiall = [phiall phi];  
    end
    
    cohpb3.cohall = cohall;
    cohpb3.phiall = phiall;
    cohpb3.f = f;

    %% coherence between binned point processes (running, binned, low-res, phase-scrambled)
    fprintf('[%d] compute coherence (running, binned, low-res, phase-scrambled)\n', sid)

    clear params cohpb4 m
    params.Fs = 50; % sampling frequency
    params.fpass = [0 25]; % frequency range of interest
    params.tapers = [5 9]; % tapers
    params.err = [2 0.05];
    params.pad = 0;
    win = 75;
    
    if size(vars.rsData1_spikect2_rnps,1)/params.Fs < win
        win = size(vars.rsData1_spikect2_rnps,1)/params.Fs;
    end
    cohall = [];
    phiall = [];
    for i = 1:nNeurons
        fprintf('processing [%d|%d]\n', i, nNeurons);
        spikect = vars.rsData1_spikect2_rnps;
        if size(spikect,1) > 10000
            spikect = spikect(2500:end-2500,:);
        end
        
        sig = spikect(:,i);
        gls = [spikect(:,1:i-1) spikect(:,i+1:end)];
        gls = nanmean(gls, 2);
    
        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr] = coherencysegpb(sig, gls, win, params, 1);
        cohall = [cohall C];
        phiall = [phiall phi];  
    end
    
    cohpb4.cohall = cohall;
    cohpb4.phiall = phiall;
    cohpb4.f = f;
    
    %% save data
    cohpb1_dat = {cohpb1.cohall, cohpb1.phiall, cohpb1.f'};
    cohpb2_dat = {cohpb2.cohall, cohpb2.phiall, cohpb2.f'};
    cohpb3_dat = {cohpb3.cohall, cohpb3.phiall, cohpb3.f'};
    cohpb4_dat = {cohpb4.cohall, cohpb4.phiall, cohpb4.f'};
    
    fpath = [sessDir, sprintf('session_%d/', sid), 'coherence_high_res.mat'];
    save(fpath, 'cohpb1_dat', 'cohpb2_dat', 'cohpb3_dat', 'cohpb4_dat');
end
