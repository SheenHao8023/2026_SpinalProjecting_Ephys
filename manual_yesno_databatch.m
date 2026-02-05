clear
close all
addpath(genpath('C:\Data_CompileAllData'));

% meta_file_dir_CIBR_
file_dir_all = {
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
% tag_name_all = {
%        'optobehavior_LSC_LSNr_CHR2_473nm_cos40Hz400msRamp_mean1.5mW_GC';...
% };

for i_session = 1:size(file_dir_all,1)  
    [suapath_root, ~, ~] = fileparts(file_dir_all{i_session});
    [session_dir, ~, ~] = fileparts(suapath_root);
    session_dir = fullfile(session_dir, '\');  
    % session_tagname = tag_name_all{i_session}; 
    
    % batch behavior .mat file
    solo_file = dir([session_dir, '*.mat']);
    if size(solo_file,1)==1 % only one file per session
        solo_filename = [session_dir, solo_file.name];
    else  % multiple solo files per session
        solo_filename = {};
        filename_tmp = dir([session_dir,'*_Session1.mat']);
        solo_filename{1} = [session_dir, filename_tmp.name];
        filename_tmp = dir([session_dir,'*_Session2.mat']);
        solo_filename{2} = [session_dir, filename_tmp.name];
    end

    % get output filename
    istr = strfind(session_dir,'\');
    animal_name = session_dir(istr(end-2)+1:istr(end-1)-1);
    date_name = [session_dir(istr(end-1)+1:istr(end)-1)];
    output_filename = [animal_name,'_',date_name];
    output_dir = session_dir(1:istr(end-1));
    
    % construct the session data object
    disp(['............ ',animal_name,'  ',date_name,' ...............']);
    % construct the behavioral data object
    obj = Solo.YesNoDiscriminationCIBRData(solo_filename); % protocol type 4/5
    % obj.Session_Name_Tag(session_tagname);
    save([output_dir,output_filename,'_behavior.mat'],'obj','-v7.3');
    
    % batch alldata .mat file
    wavesurfer_direname = dir([session_dir, '*.h5']);
    if ~isempty(wavesurfer_direname) && isscalar(wavesurfer_direname)
        wavesurfer_direname = [session_dir wavesurfer_direname.name];
        single_unit_dir1 = [session_dir, 'Kilosort_SingleUnits1\'];
        single_unit_dir2 = [session_dir, 'Kilosort_SingleUnits2\'];
        obj.load_Wavesurfer_Data(wavesurfer_direname)
        single_unit_filelist1 = dir([single_unit_dir1, 'SingleUnit*.mat']);
        single_unit_filelist2 = dir([single_unit_dir2, 'SingleUnit*.mat']);
    else
        if isempty(wavesurfer_direname)
            fprintf('未找到 .h5 文件\n');
        else
            fprintf('找到 %d 个 .h5 文件\n', length(wavesurfer_direname));
        end
        single_unit_filelist1 = [];
        single_unit_filelist2 = [];
    end
    
    % the sorted unitset 1
    Units=[];
    if size(single_unit_filelist1,1)>0
        for i_file = 1:size(single_unit_filelist1,1)
            disp('--------------------------------------------');
            disp(['processing file ',single_unit_filelist1(i_file).name]);
            load([single_unit_dir1, single_unit_filelist1(i_file).name]);
            unit.waveforms_allch=[];
            Units{i_file}=unit;
        end
    end
    obj.units{1}=Units;

    %the sorted unitset 2
    Units=[];
    if size(single_unit_filelist2,1)>0
        for i_file = 1:size(single_unit_filelist2,1)
            disp('--------------------------------------------');
            disp(['processing file ',single_unit_filelist2(i_file).name]);
            load([single_unit_dir2, single_unit_filelist2(i_file).name]);
            unit.waveforms_allch=[];
            Units{i_file}=unit;
        end
    end
    obj.units{2}=Units;

    save([output_dir,output_filename,'_allData.mat'],'obj','-v7.3');
end