%% read kilosort to singleunits
% cgs - 0 = noise / 1 = mua / 2 = good / 3 = unsorted
clear;
addpath(genpath('C:\Kilosort2-user\my_kilosort'));
Datapath={
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
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260130_1\2026-01-30_16-52-44\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260130_2\2026-01-30_18-03-25\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260131_1\2026-01-31_16-49-26\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH035\20260131_2\2026-01-31_17-27-07\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_1\2026-02-02_15-54-33\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_2\2026-02-02_16-59-51\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH037\20260202_3\2026-02-02_17-35-26\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH036\20260203_1\2026-02-03_18-11-23\';...
    'C:\Users\XinHao\Desktop\CoopXiang\CIBRXH036\20260203_2\2026-02-03_18-45-02\';...
};
% Laser_artifact={'n';...
%     'n';...
%     'n';...
%     'n';...
%     };

Probe = 2;

for k=1:length(Datapath)
    % laser_artifact=Laser_artifact{k};
    for probe=1: Probe
        for shank=1:1
            istr = strfind(Datapath{k},'\');
            datapath=[Datapath{k},'Kilosort_data',num2str(probe),num2str(shank),'\'];
            % 在session日期路径下保存single units
            [suapath_root, ~, ~] = fileparts(Datapath{k});
            [suapath_root, ~, ~] = fileparts(suapath_root);
            suapath = [suapath_root, '\Kilosort_SingleUnits', num2str(probe), '\'];  % shank体现在后文文件名中
            if ~exist(suapath, 'dir')
                mkdir(suapath);
            end
            rawdatafile=dir([datapath,'/*.bin']);
            trialtimefile=dir([datapath,'Trial_time_all*.mat']);%_new
            load([datapath,trialtimefile.name]);
            spike_sample_time = readNPY([datapath,'spike_times.npy']);
            spike_clusters = readNPY([datapath,'spike_clusters.npy']);
            spike_amplitudes = readNPY([datapath,'amplitudes.npy']);
            [cidstmp, cgs] = readClusterGroupsCSV([datapath,'cluster_group.tsv']);  %mua 1; good 2; noise 3
            [cids, channel, firing_rate]= readClusterInfoCSV([datapath,'cluster_info.tsv']);  % already +1 to be 1-384, 0-383 in kilosort phy, cids is from min to max order
            % 修改chMap使序号从1开始
            chMap = readNPY(fullfile(datapath, 'channel_map.npy')) + 1;  % Order in which data was streamed to disk; must be 1-indexed for Matlab
            chPos = readNPY(fullfile(datapath, 'channel_positions.npy')); 

            %keep only good units
            cids_good=cidstmp(find(cgs==2));
            if ~isempty(cids_good)
                % channel_good使序号从1开始
                channel_good=channel(find(ismember(cids,cids_good)));%already 1-384 range
                firing_rate_good=firing_rate(find(ismember(cids,cids_good)));
                idx=find(ismember(spike_clusters,cids_good));
                spike_clusters_good=spike_clusters(idx);
                spike_sample_time_good=spike_sample_time(idx);
                spike_amplitudes_good=spike_amplitudes(idx);
                
                gwfparams.dataDir = datapath;    % KiloSort/Phy output folder
                gwfparams.fileName = rawdatafile.name;         % .dat file containing the raw 
                gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)

                if max(chMap)<=32
                % CB64=2x32
                    fs=round(1/median(diff(Trial_time_all(1,:))));
                    gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
                    if fs==20000
                        gwfparams.wfWin = [-10 29];              % Number of samples before and after spiketime to include in waveform %[-10 29] for CB64  [-15 44] for NP 30K fs
                    else
                        gwfparams.wfWin = [-15 44];
                    end                    
                    gwfparams.nWf = 10000;                    % Number of waveforms per unit to pull out
                else
                    % NP384
                    fs=round(1/median(diff(Trial_time_all(1,:))));
                    gwfparams.nCh = 384;                     % Number of channels that were streamed to disk in .dat file
                    if fs==20000
                        gwfparams.wfWin = [-10 29];              % Number of samples before and after spiketime to include in waveform %[-10 29] for CB64  [-15 44] for NP 30K fs
                    else
                        gwfparams.wfWin = [-15 44];
                    end
                    gwfparams.nWf =1000;                    % Number of waveforms per unit to pull out
                end

                gwfparams.spikeTimes = spike_sample_time_good; % Vector of cluster spike times (in samples) same length as .spikeClusters
                gwfparams.spikeClusters = spike_clusters_good; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
                
                % 若需要修改getWaveForms内部 chMap使序号从1开始，见local func
                wf = getWaveForms(gwfparams);
                
                figure;
                for j=1:length(wf.unitIDs)
                    chj=channel_good(find(cids_good==wf.unitIDs(j)));%already 1-384 range
                    unit.ID=wf.unitIDs(j);
                    unit.firing_rate=firing_rate_good(find(cids_good==wf.unitIDs(j)));
                    chjtmp=find(chMap==chj);
                    unit.waveforms=squeeze(wf.waveForms(j,:,chjtmp,:));
                    Distance=sqrt(sum((chPos-repmat(chPos(chjtmp,:),size(chPos,1),1)).^2,2));% update by Guang 20250329
                    [ss,ii]=sort(Distance);
                    unit.waveforms_allch=squeeze(wf.waveForms(j,:,ii(1:32),:));
                    unit.wvsch=chMap(ii(1:32));
                    unit.wvspos=chPos(ii(1:32),:);

                    idx_tmp=find(spike_clusters_good==wf.unitIDs(j));
                    idx_tmp=idx_tmp(1:end-1);
                    % spike_times & trials 安全索引检查
                    spike_idx = spike_sample_time_good(idx_tmp);  % 原始索引
                    if max(spike_idx) <= size(Trial_time_all, 2)
                        unit.spike_times = Trial_time_all(1, spike_idx)';
                        unit.trials      = Trial_time_all(2, spike_idx)';
                    else
                        % 只保留合法索引
                        spike_idx = spike_idx(spike_idx <= size(Trial_time_all, 2));
                        unit.spike_times = Trial_time_all(1, spike_idx)';
                        unit.trials      = Trial_time_all(2, spike_idx)';
                    end
                    spike_idx = wf.spikeTimeKeeps(j, ~isnan(wf.spikeTimeKeeps(j,:)));
                    if max(spike_idx) < size(Trial_time_all, 2)
                        unit.spike_times_wv_keeps=Trial_time_all(1,wf.spikeTimeKeeps(j,find(~isnan(wf.spikeTimeKeeps(j,:)))));
                        unit.trials_wv_keeps=Trial_time_all(2,wf.spikeTimeKeeps(j,find(~isnan(wf.spikeTimeKeeps(j,:)))));
                        % spike索引超界问题，最后一个 spike 被截断，提前加一个条件判断
                    else
                        spike_idx = spike_idx(spike_idx <= size(Trial_time_all,2));  % 只保留合法索引
                        unit.spike_times_wv_keeps = Trial_time_all(1, spike_idx);
                        unit.trials_wv_keeps     = Trial_time_all(2, spike_idx);
                    end
                    unit.amplitudes=spike_amplitudes_good(idx_tmp);

                    if length(unit.trials)~=length(unit.amplitudes)
                       'error: spike amp mismatch';
                    end
                    if shank==1
                       unit.channel=repmat(chj,length(unit.spike_times),1);%DBC 32*2 shank1  or NP 1.0 or 2.0 %already 1-384 range
                    else
                       unit.channel=repmat(chj+32,length(unit.spike_times),1);%DBC 32*2 shank2
                    end
                    unit.stable_trials=min(unit.trials):max(unit.trials);
                    unit = func_classify_unit_intan(unit,fs,0);
                    unit.manual_quality_score='1';
                    % unit.artifact_annotation=laser_artifact;
                    if exist('oepxi_motor12')
                        unit.oepxi_motor12=oepxi_motor12;
                    else
                        unit.oepxi_motor12=[];
                    end

                    disp(['saving ',suapath, ' unit ',num2str(j)]);
                    save([suapath, 'SingleUnit',num2str(shank),'_',num2str(j),'.mat'],'unit');

                    if unit.cell_type==1
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'r');hold on;
                        subplot(2,4,2);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'r-'); hold on;
                    elseif unit.cell_type==2
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'b');hold on;
                        subplot(2,4,3);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'b-'); hold on;
                    else
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'k');hold on;
                        subplot(2,4,4);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'k-'); hold on;
                    end
                    subplot(2,4,5:8);hold on;
                    scatter(unique(unit.trials),j*ones(1,length(unique(unit.trials))),'g+');hold on;
                    xlabel('trial');
                    ylabel('unit');

                end
                saveas(gcf, [suapath, 'celltype_depth_map.png'],'png');
                close;
                clearvars -except Datapath Laser_artifact laser_artifact k position probe shank Probe
            end
        end
    end
end

%% manual correct Spike Timing kilosort
clearvars -except Datapath Probe;
close all

for probe=1:Probe
    for i_session = 1:size(Datapath,1)
        [suapath_root, ~, ~] = fileparts(Datapath{i_session});
        [suapath_root, ~, ~] = fileparts(suapath_root);
        single_unit_dir = [suapath_root,'\Kilosort_SingleUnits',num2str(probe),'\'];
        unit_files = dir([single_unit_dir,'*.mat']);
        for i_unit = 1:size(unit_files,1)
            load([single_unit_dir, unit_files(i_unit).name]); 
            unit_tmp = unit;
            isi = diff(unit_tmp.spike_times);
            i_1st_spk = diff(unit_tmp.trials)>0;
           
            i_spk_discard = find(isi<.0005  & i_1st_spk==0);

            unit_tmp.spike_times(i_spk_discard,:) = [];
            unit_tmp.amplitudes(i_spk_discard,:) = [];
            unit_tmp.trials(i_spk_discard,:) = [];
            unit_tmp.channel(i_spk_discard,:) = [];
            unit_tmp.stable_trials   =  unit.stable_trials;
            unit_tmp.cell_type = unit.cell_type;
            unit_tmp.ID=unit.ID;
            unit_tmp.firing_rate=unit.firing_rate;
            unit_tmp.manual_quality_score=unit.manual_quality_score;
            %unit_tmp.artifact_annotation=unit.artifact_annotation;

            disp([num2str(size(i_spk_discard,1)),'/',num2str(size(isi,1)+1),' spikes discarded'])
            unit = unit_tmp;
            ISI = diff(unit.spike_times);
            ISI = ISI(find(ISI<.5));
            ISI = [ISI; -ISI];
            unit.false_alarm_est = sum(abs(ISI)<.002)/length(ISI);
       
            unit.waveforms=unit.waveforms(find(~isnan(unit.waveforms(:,1))),:);
            wave_amp_tmp = range(unit.waveforms,2);
            if size(wave_amp_tmp,1)>100
                mean_wave_amp = conv(wave_amp_tmp,ones(1,100)/100,'same');
                mean_wave_amp(1:50) = mean_wave_amp(51);
                mean_wave_amp(end-49:end) = mean_wave_amp(end-50);
                wave_amp_tmp = wave_amp_tmp-mean_wave_amp;   % only look at the residues
            end
       
            mu_est = mean(wave_amp_tmp);
            sigma_est = std(wave_amp_tmp);
       
            X_min = (mu_est-sigma_est*5);
            X_max = (mu_est+sigma_est*5);
            X = X_min:(X_max-X_min)/100:X_max;
            Y_fit = normpdf(X,mu_est,sigma_est);
            Y = histc(wave_amp_tmp,X);
            Y = Y/sum(Y)/((X_max-X_min)/100);
            r = corr(Y,Y_fit');
       
            unit.miss_est = r^2;
            save([single_unit_dir, unit_files(i_unit).name],'unit'); 

        end
    end
clearvars -except Datapath session_dir position probe Probe
end

%% manual check for duplicate clusters kilosort
clearvars -except Datapath Probe;
close all

probetype=1; %1 CB64; 2 NP1.0 or NP2.0

for probe=1:Probe
    for i_session = 1:size(Datapath,1)
        [suapath_root, ~, ~] = fileparts(Datapath{i_session});
        [suapath_root, ~, ~] = fileparts(suapath_root);
        single_unit_dir = [suapath_root,'\Kilosort_SingleUnits',num2str(probe),'\'];
        if ~exist([single_unit_dir,'CheckDuplicate'])
            mkdir([single_unit_dir,'CheckDuplicate'])
        end

        unit_files = dir([single_unit_dir,'SingleUnit*.mat']);
        for i_unit = 1:size(unit_files,1)
            load([single_unit_dir, unit_files(i_unit).name]);
            unit_tmp = unit;
            unit_channel = median(unit_tmp.channel);
            unit1pos=unit.wvspos(1,:);

            for i_unit_pair = i_unit+1:size(unit_files,1)
                load([single_unit_dir, unit_files(i_unit_pair).name]);
                unit_pair_tmp = unit;
                unit_pair_channel = median(unit_pair_tmp.channel);
                unit2pos=unit.wvspos(1,:);

                if strcmp(unit_files(i_unit).name(11),unit_files(i_unit_pair).name(11))
                    if sqrt(sum((unit1pos-unit2pos).^2))<30 %um %abs(unit_channel-unit_pair_channel)<=1
                        if length(unique(unit_tmp.trials))>10&length(unique(unit_pair_tmp.trials))>10
                            close all
                            func_check_for_duplicate_unit_kilosort(unit_tmp, unit_files(i_unit).name(1:end-4), unit_pair_tmp, unit_files(i_unit_pair).name(1:end-4), single_unit_dir, probetype);
                        end
                    end
                end  
            end
        end
    end
    clearvars -except Datapath session_dir position probe probetype
end

%% local func getWaveForms
function wf = getWaveForms(gwfparams)
    % Load .dat and KiloSort/Phy output
    fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);
    filenamestruct = dir(fileName);
    dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
    nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
    wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
    mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
    chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy')) + 1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
    nChInMap = numel(chMap);
    
    unitIDs = unique(gwfparams.spikeClusters);
    numUnits = size(unitIDs,1);
    spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
    waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
    waveFormsMean = nan(numUnits,nChInMap,wfNSamples);
    xlength=size(mmf.Data.x,2);
    
    for curUnitInd=1:numUnits
        curUnitID = unitIDs(curUnitInd);
        curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
        curUnitnSpikes = size(curSpikeTimes,1);
        spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
        spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
        for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
            spikeTimekeeps_tmp=spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end);
            if spikeTimekeeps_tmp(1)>=1&&spikeTimekeeps_tmp(end)<=xlength
                tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimekeeps_tmp);
                waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
            end
        end
        waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
        disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
    end
    
    wf.unitIDs = unitIDs;
    wf.spikeTimeKeeps = spikeTimeKeeps;
    wf.waveForms = waveForms;
    wf.waveFormsMean = waveFormsMean;
end