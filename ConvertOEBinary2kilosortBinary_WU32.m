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

% Laser_trial={
%     [1 2 3];...
%     [1 2 3];...
%     [1 2 3];...
%     [1 2 3];...
%     }; %sti right/left/bilat

for k = 1:length(datapath)
    file_dir = datapath{k};
    % i_str=strfind(file_dir,'\');
    % laser_trial = Laser_trial{k};

    probe = 2; %region 1 intan headstage 1-64 ch;region 2 intan headstage 65-128 ch etc for mutiple arrays in different regions
    shank = 1; 
    mmapN=1;
    oepxi=0;
    for probe_tmp=1:probe
        for shank_tmp=1:shank
            if ~exist([file_dir,'Kilosort_data',num2str(probe_tmp),num2str(shank_tmp)]) 
                mkdir([file_dir,'Kilosort_data',num2str(probe_tmp),num2str(shank_tmp)])
            end
        end
    end
    func_process_voltage_traces_kilosort_WuElectrode(file_dir, probe, shank, mmapN);
    % master_kilosort_WuElectrode; 
end

%% Local func: func_process_voltage_traces_kilosort_WuElectrode
function func_process_voltage_traces_kilosort_WuElectrode(file_dir, probe, shank, mmapN)
    oebin_file = dir(fullfile(file_dir, '**', 'structure.oebin'));
    oebin_path = fullfile(oebin_file(1).folder, oebin_file(1).name);
    D = load_open_ephys_binary(oebin_path, 'continuous', mmapN, 'mmap');
    fs=D.Header.sample_rate;

    probe1_channels = arrayfun(@(x) ['CH' num2str(x)], 17:48, 'UniformOutput', false);
    probe2_channels = arrayfun(@(x) ['CH' num2str(x)], 81:112, 'UniformOutput', false);
    channel_index = cell(1, probe);   

    if probe == 1
        if D.Header.num_channels == 72 % 如果记录了全部通道，则删除17~48以外的无效通道
            delete_list = [ ...
                arrayfun(@(x) sprintf('CH%d', x), 1:16,  'UniformOutput', false), ...
                arrayfun(@(x) sprintf('CH%d', x), 49:64, 'UniformOutput', false) ...
            ];
            idx_del = ismember({D.Header.channels.channel_name}, delete_list);
            D.Header.channels(idx_del) = [];
            D.Header.num_channels = numel(D.Header.channels);
            D.mapped_copy = int16(D.Data.Data(1).mapped);
            D.mapped_copy(idx_del, :) = [];
        else
            D.mapped_copy = D.Data.Data(1).mapped;
        end
        channel_idx{1} = find(ismember({D.Header.channels.channel_name}, probe1_channels)); % number vector, 如果少记录通道
    else  % probe = 2
        if D.Header.num_channels == 136 % 如果记录了全部通道，则删除17~48, 81~112以外的无效通道
            delete_list = [ ...
                arrayfun(@(x) sprintf('CH%d', x), 1:16,  'UniformOutput', false), ...
                arrayfun(@(x) sprintf('CH%d', x), 49:80, 'UniformOutput', false) ...
                arrayfun(@(x) sprintf('CH%d', x), 113:128, 'UniformOutput', false) ...
            ];
            idx_del = ismember({D.Header.channels.channel_name}, delete_list);
            D.Header.channels(idx_del) = [];
            D.Header.num_channels = numel(D.Header.channels);
            D.mapped_copy = int16(D.Data.Data(1).mapped);
            D.mapped_copy(idx_del, :) = [];
        else
            D.mapped_copy = D.Data.Data(1).mapped;
        end
        channel_idx{1} = find(ismember({D.Header.channels.channel_name}, probe1_channels));
        channel_idx{2} = find(ismember({D.Header.channels.channel_name}, probe2_channels)) - numel(channel_idx{1}); % number vector, 如果少记录通道
    end

    probe_ch = (D.Header.num_channels-8) / probe;
    shank_ch = 32;

    wavesurfer_data= D.mapped_copy(D.Header.num_channels-7,:); % BNC1 ch65
    bitcode_data= D.mapped_copy(D.Header.num_channels-6,:); % BNC2 ch66
    Ephystrig_data= D.mapped_copy(D.Header.num_channels-5,:); % BNC3 ch67, intan trigger
    LickingLeft_data = D.mapped_copy(D.Header.num_channels-2,:); % BNC6
    LickingRight_data = D.mapped_copy(D.Header.num_channels-1,:); % BNC7
    wavesurfer_data(1:30000)=0;
    bitcode_data(1:30000)=0;
    Ephystrig_data(1:30000)=0;
    LickingLeft_data(1:30000)=0;
    LickingRight_data(1:30000)=0;

    wavesurfer_data=wavesurfer_data*D.Header.channels(end).bit_volts;
    bitcode_data=bitcode_data*D.Header.channels(end).bit_volts;
    Ephystrig_data=Ephystrig_data*D.Header.channels(end).bit_volts;
    LickingRight_data = LickingRight_data*D.Header.channels(end).bit_volts;
    LickingLeft_data = LickingLeft_data*D.Header.channels(end).bit_volts;

    trial_on_tmp= find(diff(Ephystrig_data>1.5)==1)+1;
    trial_on(1)=trial_on_tmp(1);
    trial_on=[trial_on(1) trial_on_tmp(find(diff(trial_on_tmp)>2*fs)+1)];
    trial_off_tmp= find(diff(Ephystrig_data>1.5)==-1)+1;
    trial_off(1)=trial_off_tmp(1);
    trial_off=[trial_off(1) trial_off_tmp(find(diff(trial_off_tmp)>2*fs)+1)];

    for probe_tmp=1:probe
        for shank_tmp=1:shank
            TrialNum_prev = 0;
            TrialNum = 0;
            Trial_time_all=[];
            file_dir_output=[file_dir,'kilosort_data',num2str(probe_tmp),num2str(shank_tmp),'\'];
            output_file_kilosort= [file_dir_output,'kilosort_probe',num2str(probe_tmp),'shank',num2str(shank_tmp),'_data.bin'];
            fidout = fopen(output_file_kilosort, 'w'); 
            for trial_k=1:length(trial_on)
                if trial_k==1
                    sample_point=1:trial_off(trial_k);       
                elseif trial_k>1&& trial_k<length(trial_on)
                    sample_point=trial_off(trial_k-1)+1:trial_off(trial_k);
                else
                    sample_point=trial_off(trial_k-1)+1:length(wavesurfer_data);
                    if length(sample_point)/fs>301 
                        sample_point=trial_off(trial_k-1)+1:trial_off(trial_k-1)+300*fs; 
                    end
                end
                TimeStamps=(sample_point-trial_on(trial_k))/fs;
                bitcode_trial=bitcode_data(sample_point);

                TrialNum_prev = TrialNum;
                try
                    TrialNum = func_read_bitcode(bitcode_trial,TimeStamps);
                    disp(['done!  Matched to solo trial #',num2str(TrialNum)]);
                catch
                    TrialNum = TrialNum_prev+1;
                    warning(['Bitcode Failed!  Assigned to solo trial #',num2str(TrialNum)]);
                end

                Trial_time_tmp=[TimeStamps;ones(1,length(TimeStamps))*TrialNum];
                Trial_time_all=[Trial_time_all Trial_time_tmp];
    
                ch_MUA = [];
                ch_MUA = D.mapped_copy((probe_tmp-1)*probe_ch+1:probe_tmp*probe_ch,sample_point)*D.Header.channels(1).bit_volts;
                ch_MUA = double(ch_MUA');
                ch_MUA = ch_MUA(:,(shank_tmp-1)*shank_ch+1:shank_tmp*shank_ch); 
                for i_ch = 1:length(channel_idx{probe_tmp}) % bandpass filter
                    ch_tmp = timeseries(ch_MUA(:, channel_idx{probe_tmp}(i_ch)),TimeStamps);
                    ch_tmp_MUA = idealfilter(ch_tmp, [300 6000], 'pass');
                    ch_MUA(:, channel_idx{probe_tmp}(i_ch)) = ch_tmp_MUA.data;
                end

                % substract off common noise
                commonNoise = trimmean(ch_MUA,40,2);
                i_noise=[find(commonNoise>150);find(commonNoise<-200)]; 
                t_sample=1.57;
                for i_ch=1:length(channel_idx{probe_tmp}) % change the time here acorrding to your behavior and photostim
                    t_post_stim = [0 t_sample+1.3]; 
                    i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
                    idx1=i_post_stim(end);
                    X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
                    b1 = regress(ch_MUA(i_post_stim,i_ch),X);
                    t_post_stim = [t_sample+1.3 t_sample+2.1]; 
                    i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
                    idx2=i_post_stim(end);
                    X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
                    b2 = regress(ch_MUA(i_post_stim,i_ch),X);
                    t_post_stim = [t_sample+2.1 t_sample+4.5]; 
                    i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
                    X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
                    b3 = regress(ch_MUA(i_post_stim,i_ch),X);
                    ch_MUA(:,i_ch) = ch_MUA(:,i_ch) - [commonNoise(1:idx1)*b1(2);commonNoise(idx1+1:idx2)*b2(2);commonNoise(idx2+1:end)*b3(2)];
                end
                if ~isempty(i_noise)
                    for j=1:length(i_noise)
                        if i_noise(j)>1 && i_noise(j)<length(commonNoise)-1
                            ch_MUA(i_noise(j)-1:i_noise(j)+1,:)=zeros(3,size(ch_MUA,2));
                        end
                    end
                end
                clear commonNoise

                % if ismember(TrialNum,laser_trial)
                %      output_file_name = [file_dir_output,'raw_trace_probe',num2str(probe_tmp),'shank',num2str(shank_tmp),'_laser_trial_',num2str(TrialNum),'.mat'];
                %      disp(['saving laser trial: ',output_file_name]);
                %      save(output_file_name,'ch_MUA','TimeStamps','-v7.3');
                % end
                % disp(['write: trial ',num2str(TrialNum)]);
                ch_MUA = int16(ch_MUA');
                fwrite(fidout, ch_MUA, 'int16');
                clear ch_MUA ch_tmp ch_tmp_MUA TimeStamps Trial_time_tmp

            end
            fclose(fidout);
            output_file=[file_dir_output,'Trial_time_all_probe',num2str(probe_tmp),'shank',num2str(shank_tmp),'_trial_',num2str(Trial_time_all(2,1)),'_',num2str(Trial_time_all(2,end)),'.mat'];
            save(output_file,'Trial_time_all','-v7.3');
            events=find(Trial_time_all(1,:)==0)'/fs;
            writematrix(events, fullfile(file_dir_output, 'events.csv'));
        end
    end
return
end