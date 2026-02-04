%% mannually put in the duplicate unit pair (merge single unit for each subfolder if have)
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
Probe = 2;

single_unit_dir='C:\Users\XinHao\Desktop\DiffLearn_MultiRegion\Ephys\CIBRXH027\2025_11_25\2025-11-25_12-29-37\Kilosort_SingleUnits1\';
combine_pairs={
   };

if size(combine_pairs,1)>0    

    for i_pair = 1:size(combine_pairs,1)
        tf=0;
        pair1=combine_pairs{i_pair,1};
        pair2=combine_pairs{i_pair,2};
        if exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tf=1;
        elseif ~exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tmp=[];
            for k=1:size(combine_pairs,1)
                tmp(k)=strcmp(combine_pairs{k,2},pair1);
            end
            tmp=find(tmp);
            pair1=combine_pairs{tmp(1),1};
            tf=1;
        elseif exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&~exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tmp=[];
            for k=1:size(combine_pairs,1)
                tmp(k)=strcmp(combine_pairs{k,2},pair2);
            end
            tmp=find(tmp);
           pair1=combine_pairs{tmp(1),1};
           tf=1;
        end
        if tf
            disp(['combining SingleUnit ',pair1,'  ',pair2]);
            load([single_unit_dir, 'SingleUnit',pair1,'.mat']);
            unit1 = unit;
            load([single_unit_dir, 'SingleUnit',pair2,'.mat']);
            unit2 = unit;
            clear unit;
            [unit] = func_combine_duplicate_unit_kilosort(unit1,unit2);
            isi = diff(unit.spike_times);
            i_1st_spk = diff(unit.trials)>0;
            i_spk_discard = find(isi<.0005  & i_1st_spk==0);
            
            unit.spike_times(i_spk_discard,:) = [];
            unit.amplitudes(i_spk_discard,:) = [];
            unit.trials(i_spk_discard,:) = [];
            unit.channel(i_spk_discard,:) = [];
            
            disp([num2str(size(i_spk_discard,1)),'/',num2str(size(isi,1)+1),' spikes discarded'])
            
            % score unit quality
            ISI = diff(unit.spike_times);
            ISI = ISI(ISI<.5);
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
            
            disp(['Updateing SingleUnit', pair1]);
            save([single_unit_dir, 'SingleUnit',pair1,'.mat'], 'unit');
    
            disp(['Deleting SingleUnit', pair2]);
            delete([single_unit_dir, 'SingleUnit',pair2,'.mat']);
            if exist([single_unit_dir, 'SingleUnit',pair2,'.png'],'file')
               delete([single_unit_dir, 'SingleUnit',pair2,'.png']);
            end
        end
    end
end

%% manual final plot and check badunits kilosort
clearvars -except Datapath Probe;
close all

for probe=1:Probe
    for i_session = 1:size(Datapath,1)
        [suapath_root, ~, ~] = fileparts(Datapath{i_session});
        [suapath_root, ~, ~] = fileparts(suapath_root);
        single_unit_dir = [suapath_root,'\Kilosort_SingleUnits',num2str(probe),'\'];
        if exist(single_unit_dir)
            if ~exist([single_unit_dir,'BadUnits'])
                mkdir([single_unit_dir,'BadUnits'])
            end
            unit_files = dir([single_unit_dir,'*.mat']);
    
            for i_unit = 1:size(unit_files,1)
                load([single_unit_dir, unit_files(i_unit).name]); 
                figure;
                subplot(2,3,1);
                isi = diff(unit.spike_times);
                isi = isi(find(isi<.5));
                isi = [isi; -isi];
                edges = -.03:.00025:.03;
                n = histcounts(isi, edges);       % 长度为 length(edges)-1
                n(end+1) = sum(isi == edges(end)); 
                plot(edges, n, 'r')
                if max(n)~=0
                    axis([-.02 .02 0 max(n)]);
                end
                xlabel('s');
                ylabel('count');
                subplot(2,3,4);
                if unit.cell_type == 1
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'r')
                end
                if unit.cell_type == 2
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'b')
                end
                if unit.cell_type == 0
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'k')
                end
                xlabel('sample count');
                ylabel('uV');
       
                spike_times_psth = {};
                n_trial = 0;
                trial0=min(unit.trials):max(unit.trials);
                trial_tmp=intersect(1:999,trial0)'; % only behavior trial

                if size(trial_tmp,1)>1
                    trial_tmp=trial_tmp';
                end
                if ~isempty(trial_tmp)
                    for i_trial = trial_tmp
                        n_trial = n_trial+1;
                        spike_times_psth{n_trial,1} = unit.spike_times(unit.trials==i_trial)';
                    end
                end
                if size(spike_times_psth,1)>10&max(unit.spike_times)>0
                    subplot(4,3,5);
                    [psth, t] = func_getPSTH(spike_times_psth,0,max(unit.spike_times));
                    bar(t,psth,'k');hold on;
                    xlim([0 min([10 max(unit.spike_times)])])
                    xlabel('s');
                    ylabel('spikes/s');
                    subplot(4,3,2);
                    trial_idx=find(ismember(unit.trials,trial_tmp));
                    plot(unit.spike_times(trial_idx),unit.trials(trial_idx),'.k');
                    xlim([0 min([10 max(unit.spike_times)])])
                    xlabel('s');
                    ylabel('trial');
                end
                subplot(4,3,[3 6]);
                trial_idx=find(unit.trials>=1000);
                plot(unit.spike_times(trial_idx),unit.trials(trial_idx),'.k');hold on;
                xlim([-10 10])
                xlabel('s');
                ylabel('trial');

                saveas(gcf,[single_unit_dir, unit_files(i_unit).name(1:end-4),'.png'],'png');
                close;
            end
        end
    end
end

%% remove bad bitcode trial
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
           idx=unique(unit.trials);
           idx=find(unit.trials<=idx(5)|unit.trials>=idx(end-4));
           unit.spike_times(idx)=[];
           unit.trials(idx)=[];
           unit.amplitudes(idx)=[];
           unit.channel(idx)=[];
           unit.stable_trials=min(unit.trials):max(unit.trials);
           save([single_unit_dir, unit_files(i_unit).name],'unit');
        end
    end
end
