clear
addpath(genpath('C:\Kilosort2-user\my_kilosort'));
datapath={
    % CIBRXH027, probe = 1
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_11_25\2025-11-25_12-29-37\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_11_26\2025-11-26_10-05-43\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_11_27\2025-11-27_10-06-00\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_11_28\2025-11-28_12-27-55\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_12_06\2025-12-06_09-48-52\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_12_07\2025-12-07_09-51-21\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_12_08\2025-12-08_14-16-55\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_12_09\2025-12-09_09-47-14\';... 
    % 'C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_12_11\2025-12-11_10-45-22\';... 

    % CIBRXH035/36/37, probe = 2
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260130_1\2026-01-30_16-52-44\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260130_2\2026-01-30_18-03-25\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260131_1\2026-01-31_16-49-26\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260131_2\2026-01-31_17-27-07\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_1\2026-02-02_15-54-33\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_2\2026-02-02_16-59-51\';...
    % 'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_3\2026-02-02_17-35-26\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH036\20260203_1\2026-02-03_18-11-23\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH036\20260203_2\2026-02-03_18-45-02\';...
};

for k = 1:length(datapath)
    probe = 2; 
    shank = 1; 
    for probe_tmp=1:probe
        for shank_tmp=1:shank
            rootZ = [datapath{k},'Kilosort_data',num2str(probe_tmp),num2str(shank_tmp),'\']; %,num2str(position) the raw data binary file is in this folder
            addpath(genpath('C:\Kilosort2-master')) % path to kilosort folder
            addpath('C:\Kilosort2-user\my_kilosort\npy-matlab-master\npy-matlab') % for converting to Phy
            rootH = rootZ; % path to temporary binary file (same size as data, should be on fast SSD)
            chanMapFile = 'C:\Users\XinHao\Desktop\MDL\Ephys\WU32_kilosortChanMap.mat';
            % 即使OE记录通道未从1开始，chanmap构建必须从1开始编号

            ops.trange = [0 Inf]; % time range to sort
            ops.NchanTOT    = 32; % total number of channels in your recording for each probe
            
            % pathToYourConfigFile = 'C:\Kilosort2-user\my_kilosort\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
            % run(fullfile(pathToYourConfigFile, ['configFile32',num2str(shank_ch),'.m']))
            ops.chanMap = fullfile(chanMapFile);
            % ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file
            trialtimefile=dir([rootZ,'Trial_time_all*.mat']); 
            load([rootZ,trialtimefile.name]);
            samplerate=round(1/median(diff(Trial_time_all(1,:)))); % 自动识别采样率
            ops.fs = samplerate;  
            ops.fshigh = 300; % frequency for high pass filtering
            ops.minfr_goodchannels = 0.1;  % minimum firing rate on a "good" channel (0 to skip)
            ops.Th = [10 4];  % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
            ops.lam = 10;  % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
            ops.AUCsplit = 0.9;  % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
            ops.minFR = 1/50;  % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
            ops.momentum = [20 400];  % number of samples to average over (annealed from first to second value) 
            ops.sigmaMask = 30; % spatial constant in um for computing residual variance of spike
            ops.ThPre = 8; % threshold crossings for pre-clustering (in PCA projection space)
            % danger, changing these settings can lead to fatal errors
            % options for determining PCs
            ops.spkTh           = -4;%-6;      % spike threshold in standard deviations (-6)
            ops.reorder         = 1;       % whether to reorder batches for drift correction. 
            ops.nskip           = 25;  % how many batches to skip for determining spike PCs
            ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
            % ops.Nfilt               = 1024; % max number of clusters
            ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
            ops.whiteningRange      = 32; % number of channels to use for whitening each channel
            ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.nPCs                = 3; % how many PCs to project the spikes into
            ops.useRAM              = 0; % not yet available
            ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    
            % channel map file check
            fs = dir(fullfile(rootZ, 'chan*.mat'));
            if ~isempty(fs)
                ops.chanMap = fullfile(rootZ, fs(1).name);
            end
    
            % find the binary file
            fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
            ops.fbinary = fullfile(rootZ, fs(1).name);
    
            % preprocess data to create temp_wh.dat
            rez = preprocessDataSub(ops);
    
            % time-reordering as a function of drift
            rez = clusterSingleBatches(rez);
    
            % saving here is a good idea, because the rest can be resumed after loading rez
            save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');
    
            % main tracking and template matching algorithm
            rez = learnAndSolve8b(rez);
    
            % final merges
            rez = find_merges(rez, 1);
    
            % final splits by SVD
            rez = splitAllClusters(rez, 1);
    
            % final splits by amplitudes
            rez = splitAllClusters(rez, 0);
    
            % decide on cutoff
            rez = set_cutoff(rez);
            fprintf('found %d good units \n', sum(rez.good>0))
    
            % write to Phy
            fprintf('Saving results to Phy  \n')
            rezToPhy(rez, rootZ);
    
            clearvars -except datapath Laser_trial k position probe shank probe_tmp shank_tmp mmapN
        end
    end
end