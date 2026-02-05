% XH 2025/12/28 训练全阶段
% ZG, 3/12.
% GuangChen 10/8/23


classdef YesNoDiscriminationCIBRData < handle
    
    properties
        %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%%
        % trim = [10 10]; % The number of trials to trim from beginning and end.
        trim; % The number of trials to trim from beginning and end.
        
        mouseName = '';
        sessionName = '';
        sessionType = '';
        sessionNameTag = '';
        
        eventsHistory;  % RTLSM state history
        stimType;       % sample, delay or mixted et al
        stimProb;       % 0.25 0.5 0.75 or any other number
        sides;          % right or left
        sidesTypes;     % s: stim; n: no-stim
        PoleSideType;   % Pole right or left for bilateral task
        contextType;
        tasteType;

        TrialPoleRight;
        TrialPoleLeft;
        TrialSideTypes;

        photo_input1History;% in unit V, can be translated into real laser power
        photo_input2History;  % in unit V, can be translated into real coordinates
        photo_input3History; % input 1 2 3 depend on experiments, could be power, x y coordinates or opto switch channel
        stimStateHistory;
        
        delayDuration;
      
        trials;         % structure contains basic statistics of that trial type, trials.hitHistory, trials.missHistory, trials.noResponseHistory
        TrialStartTimestamp;
        TrialEndTimestamp;

        photo_input1;    % 4xnx2 matrix; right: 4xnx1; left: 4xnx2; (1,:,:) hit rate; (2,:,:) correct number; (3,:,:) totoal number; (4,:,:) aom voltage;        
        photo_input2;
        photo_input3;
        video_trigger_input1;
        video_trigger_input2;
        video_trigger_input3;
        video_trigger_input4;
        
        licking;    % structure has licking.response_time, licking.time, licking.side
        poleTime; % n X 2 vector, state 40 to be 0, use state 41 58 for pole moving time
        cueTime;   % state 56
        
        
        %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%%
        soloTrialIndex; % This variable only becomes relevant when the arrays were built from multiple solo files
                        % i.e. sessions w/ 2 solo files
                        % [col1: which file the trial were from, e.g. 1 or 2 to indicate file1 or file2
                        %  col2: what is the trial number of this trial in its original solo file]

        
        wavesurfer;   % laser trace, lick trace, bit code, frame trigger, frame gate
        wavesurfer_fs   % Wavesurfer sampling rate
        
        units;
        user_var
        unit_xcoords;
        unit_ycoords;
    end
    
%     properties (Dependent = true)
%         
%         
%         multiGoPosition;
%         positionRange;
%         
%         seqFileOffset;
%         
%         trimmedTrialNums;
%         volumeH20Consumed;
%         EmbCLog;
%         stimGiven;
%         numStims;
%         
%         contact
%     end
    
    methods (Access = public)
        function obj = YesNoDiscriminationCIBRData(fileName, trimData, wavesurfer_fs)
            %
            % function obj = BTStimTimingArray(x, session_name)
            %
            % Input argument 'x' is either BPOD file name string
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin >= 3 && ~isempty(wavesurfer_fs)
                obj.wavesurfer_fs = wavesurfer_fs;
            end
            
            %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%%
            trimData = [10 10]; % The number of trials to trim from beginning and end. %GC change trimData to trim
            if iscell(fileName)
                
                for i_solo_file = 1:length(fileName)
                    if i_solo_file == 1
                        obj = Solo.YesNoDiscriminationCIBRData(fileName{i_solo_file}, [], obj.wavesurfer_fs);
                    else
                        obj_tmp = Solo.YesNoDiscriminationCIBRData(fileName{i_solo_file}, [], obj.wavesurfer_fs);
                        obj.addYesNoDiscriminationCIBRData(obj_tmp);
                        clear obj_tmp
                    end
                end
                
                
                
            elseif isstr(fileName)
                %%%   read raw solo file
                load(fileName);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  basic information: session type, mouse name, file name, mouse weight
                i_str = findstr(fileName,'_');
                obj.mouseName = fileName(1:(i_str(1)-1));
                obj.sessionName = fileName;
                if isfield(SessionData,'TrialStartTimestamp')
                obj.TrialStartTimestamp=SessionData.TrialStartTimestamp';
                obj.TrialEndTimestamp=SessionData.TrialEndTimestamp';
                else
                    obj.TrialStartTimestamp=nan(SessionData.nTrials,1);
                    obj.TrialEndTimestamp=nan(SessionData.nTrials,1);
                end

                tmp = [SessionData.TrialSettings.GUI];
                tmp = [tmp.ProtocolType];
                sessionType = {};
                if isfield(SessionData,'Settings')
                    for i=1:SessionData.nTrials
                    sessionType{i,1} = SessionData.Settings.GUIMeta.ProtocolType.String{tmp(i)};
                    end
                else
                for i=1:SessionData.nTrials
                    sessionType{i,1} = SessionData.SettingsFile.GUIMeta.ProtocolType.String{tmp(i)};
                end
                end
                obj.sessionType=sessionType;       % NL 8/5/13
                clear tmp i sessionType
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  history: state machine states, stim prob, stim type, trial sides, trial types, AOM power, xGalvo, yGalvo, stim state
                
                obj.eventsHistory=(SessionData.RawEvents.Trial)';
                
                % stimulation flag
                try
                    tmp = [SessionData.TrialSettings.GUI];
                    obj.stimProb = [tmp.StimProbe]';
                    obj.stimType=SessionData.StimTypes; 
                catch
                    obj.stimProb = [];
                    obj.stimType = [];
                end
                
                %%-----need to add---->obj.stimType=saved_history.TrialTypeSection_StimType;
                if isfield(SessionData,'TrialTypes')
                    trial_type_tmp = SessionData.TrialTypes;      %0's (right) or 1's (left)
                    trial_type_tmp(trial_type_tmp==0) = 'r';
                    trial_type_tmp(trial_type_tmp==1) = 'l';
                    obj.sides = trial_type_tmp';
                    clear trial_type_tmp;
                end
                if isfield(SessionData,'TrialSideTypes')
                trial_type_tmp = SessionData.TrialSideTypes;      %0's (pole right) or 1's (pole left)
                trial_type_tmp(trial_type_tmp==0) = 'R';
                trial_type_tmp(trial_type_tmp==1) = 'L';
                obj.PoleSideType = trial_type_tmp';
                clear trial_type_tmp;
                end

                % for attention
                if isfield(SessionData,'TrialPoleRight')
                    obj.TrialPoleRight = SessionData.TrialPoleRight'; %0's (posterior) or 1's (anterior)
                end
                if isfield(SessionData,'TrialPoleLeft')
                    obj.TrialPoleLeft = SessionData.TrialPoleLeft'; %0's (posterior) or 1's (anterior)
                end

                
                if isfield(SessionData,'ContextTypes')
                    obj.contextType=SessionData.ContextTypes';
                end
                if isfield(SessionData,'TasteTypes')
                    obj.tasteType=SessionData.TasteTypes';
                end

                %%-----need to add---->obj.types = saved.TrialTypeSection_previous_types';
                
                % photostimulation type
                try
                    sidesTypes_tmp = {};
                    for x = 1:SessionData.nTrials
                        if SessionData.StimTypes(x)==1
                            sidesTypes_tmp{x,1} = 'sample_period';
                        elseif SessionData.StimTypes(x)==2
                            sidesTypes_tmp{x,1} = 'delay_period';
                        elseif SessionData.StimTypes(x)==3
                            sidesTypes_tmp{x,1} = 'response_period';
                        elseif SessionData.StimTypes(x)==4
                            sidesTypes_tmp{x,1} = 'contextbaseline_trial';
                        elseif SessionData.StimTypes(x)==0
                            sidesTypes_tmp{x,1} = '_n';
                        end
                    end
                    obj.sidesTypes =  sidesTypes_tmp;
                    clear sidesTypes_tmp;
                catch
                    obj.sidesTypes = [];
                end
                    
                
                obj.photo_input1History=[]; %<-------------- need to fill in from Ephys
                %%-----need to add---->obj.xGalvoHistory=saved.stim_timing2AFCobj_stim_xGalvo_pos_history;
                %%-----need to add---->obj.yGalvoHistory=saved.stim_timing2AFCobj_stim_yGalvo_pos_history;
                %%-----need to add---->obj.Stim_State_History;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%% history: hit, miss, no response for total trials, right trials, left trials, right stim trials, left stim trials
                Outcomes = []; EarlyLicks = [];
                for x = 1:SessionData.nTrials
                    if SessionData.TrialSettings(x).GUI.ProtocolType==3 || SessionData.TrialSettings(x).GUI.ProtocolType==4 || SessionData.TrialSettings(x).GUI.ProtocolType==5 & isempty(findstr(obj.sessionName,'passive'))
                        if ~isnan(SessionData.RawEvents.Trial{x}.States.Reward(1))
                            Outcomes(x) = 1;    % correct
                        elseif ~isnan(SessionData.RawEvents.Trial{x}.States.TimeOut(1))
                            Outcomes(x) = 0;    % error
                        elseif ~isnan(SessionData.RawEvents.Trial{x}.States.NoResponse(1))
                            Outcomes(x) = 2;    % no repsonse
                        else
                            Outcomes(x) = 3;    % others
                        end
                    else
                        Outcomes(x) = 3;        % others
                    end
                    
                    if SessionData.TrialSettings(x).GUI.ProtocolType==5 & isempty(findstr(obj.sessionName,'passive'))
                        if ~isnan(SessionData.RawEvents.Trial{x}.States.EarlyLickSample(1)) | ~isnan(SessionData.RawEvents.Trial{x}.States.EarlyLickDelay(1))
                            EarlyLicks(x) = 1;    % lick early
                        else
                            EarlyLicks(x) = 0;    % others
                        end
                    else
                        EarlyLicks(x) = 0;        % others
                    end
                end
                
                obj.trials.hitHistory=(Outcomes==1)';
                obj.trials.missHistory=(Outcomes==0)';
                obj.trials.noResponseHistory=(Outcomes==2)';
                obj.trials.EarlyLicksHistory=(EarlyLicks==1)';
                obj.trials.trialNums = SessionData.nTrials;
                
                
                clear x Outcomes EarlyLicks
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%% trim history: trim the beginning, the end trials or any given trials
                
%                 if(exist('trimData'))
                    if(length(trimData)==2)
                        trimTrials=[1:trimData(1) SessionData.nTrials-trimData(2)+1:SessionData.nTrials];
                    elseif(length(trimData)>2)
                        trimTrials=trimData;
                    end
%                 else
%                     % choose only the stim/performing periods
%                     StimTrials_tmp = (obj.stimProb~=0);
%                     i_performing = find(StimTrials_tmp==1);
%                     if ~isempty(i_performing)
%                         seg_break_pt = i_performing(diff(i_performing)>1);
%                         seg_break_pt = [seg_break_pt; i_performing(end)];
%                         for i_tmp = seg_break_pt'
%                             if i_tmp<6
%                                 StimTrials_tmp(1:i_tmp) = 0;
%                             else
%                                 StimTrials_tmp(i_tmp-5:i_tmp) = 0;
%                             end
%                         end
%                     end
%                     trimTrials = (StimTrials_tmp==0);
                    
%                 end
                obj.trim = trimTrials';
                
                % trims trials (hit miss no response), trialsRNoStim, trialsLNoStim, trialsRStim, trialsLStim, aomPowerHistory, xGalvoHistory, yGalvoHistory
%                 obj.Trim_Trials(trimTrials);
                                
                obj.soloTrialIndex = [ones(SessionData.nTrials,1) (1:SessionData.nTrials)'];
                
            else
                error('input data type not recognized');
            end
            
            
        end
        
        function Session_Name_Tag( obj, nameTag )
            % nameTag is a string with contains a lot of attributes describing the session, here are some examples. You can permuate the key words to
            % account for complex situations
            %             % default stimualtion protocol, 'cosine', if no 'pulse' or 'continuous'
            %
            %             TASK_AREA_STIMCONDITION
            %
            %             TASK:
            %
            %             TwoPoleNoDelay
            %             TwoPoleWithDelay
            %             MultiPoleNoDelay
            %             ReversePoleNoDelay
            %             ReversePoleWithDelay
            %             TwoPoleRecordingWithDelay
            %             PassiveRecording
            %
            %             AREA:
            %
            %             S1
            %             leftALM
            %             rightALM
            %             MappingS1
            %             MappingWholeBrain
            %             bilateralALM
            %
            %             STIMCONDITION
            %
            %             singlePower_SampleDelay
            %             singlePower_Sample
            %             singlePower_Delay
            %             singlePower_Reward
            %             multiplePower_SampleDelay
            %             multiplePower_Sample
            %             multiplePower_Delay
            %             multiplePower_Reward            
            %             
            %             STIMPROTOCOL
            %             cosine
            %             continuous
            %             pulse
            %             protocolComparison for passive recording and behavior sessions comparing effects of silencing
            %             spatialMapping for passive recoring in mapping the spatial extent of stimulation
            %
            %             RESEARCHER
    	    %             research last name, to facilitate comparison of data from different people
            %
            %             e.g.:
            %             TwoPoleWithDelay_MappingWholeBrain_singlePower_Sample  -- whole brain mapping session
            %             TwoPoleWithDelay_leftALM_multiplePower_SampleDelay  -- left ALM dose response curve for sample and delay
            
            obj.sessionNameTag=nameTag;
            
        end
        
        function Alarm_Nums(obj) % change state alarm, alarm history
            
            %             which_trials='';         obj.Alarm_Which_Trials(which_trials);
            %             which_trials='RNoStim';        obj.Alarm_Which_Trials(which_trials);
            %             which_trials='LNoStim';        obj.Alarm_Which_Trials(which_trials);
            %             which_trials='RStim';    obj.Alarm_Which_Trials(which_trials);
            %             which_trials='LStim';    obj.Alarm_Which_Trials(which_trials);
            
            
            obj.trials.alarmNums = find(obj.trials.EarlyLicksHistory(1:obj.trials.trialNums(1),1)==1);


        end
                
        function Pole_Time(obj)

            obj.poleTime=zeros(sum(obj.trials.trialNums(1)),2);
            for ntrial=1:sum(obj.trials.trialNums(1))
                State_Matrix=obj.eventsHistory{ntrial}.States;
                
                obj.poleTime(ntrial,:)=[State_Matrix.SamplePeriod(1,1) State_Matrix.SamplePeriod(end,2)];
            end
            
        end
        
        function Cue_Time(obj)
            
            
            obj.cueTime=zeros(sum(obj.trials.trialNums(1)),2);
            for ntrial=1:sum(obj.trials.trialNums(1))
                State_Matrix=obj.eventsHistory{ntrial}.States;
                
                obj.cueTime(ntrial,:)=[State_Matrix.ResponseCue(1,1)];
            end
            
        end
        
        
        function  [Licking_Side, Licking_Timing, Response_Time]=Lick_SideTime( obj );
            
            %   the function identifies licking events from every trial during a behavior session
            %
            %   Output: 
            %           Licking_Side, indicate which side the animal first licks 
            %           Licking_Timing, [lick_time lick_side]
            %           Response_Time,  from onset of cue to first lick
            %
            % 10/3/2011 by Zengcai Guo
            
            if( length(obj.eventsHistory) ==0)
                disp( 'must present a state matrix trial events' )
                return;
            end
            
            Licking_Side = [];
            Licking_Timing = {};
            Response_Time = [];
            for i_trial = 1:sum(obj.trials.trialNums(1))

                eventHistory_iTrial = obj.eventsHistory{i_trial};
                
                lick_data_tmp = [];
                if isfield(eventHistory_iTrial.Events,'Port1In')
                    lick_data_tmp = [(eventHistory_iTrial.Events.Port1In)' (eventHistory_iTrial.Events.Port1In)'*0+'l'];
                end
                if isfield(eventHistory_iTrial.Events,'Port2In')
                    lick_data_tmp = [lick_data_tmp;[(eventHistory_iTrial.Events.Port2In)' (eventHistory_iTrial.Events.Port2In)'*0+'r']];
                end
                if ~isempty(lick_data_tmp)
                lick_data_tmp = sortrows(lick_data_tmp);
                Licking_Timing{i_trial,1} = lick_data_tmp;
                else
                    Licking_Timing{i_trial,1} = nan;
                end
                
                try                    
                    cue_time_tmp = [eventHistory_iTrial.States.ResponseCue(1,1)];
                    if ~isempty(lick_data_tmp)
                        lick_data_tmp = lick_data_tmp(lick_data_tmp(:,1)>cue_time_tmp,:);
                        Licking_Side(i_trial,1) = lick_data_tmp(1,2);
                        Response_Time(i_trial,1) = lick_data_tmp(1,1)-cue_time_tmp;
                    else
                        Licking_Side(i_trial,1) = nan;
                        Response_Time(i_trial,1) = nan;
                    end
                catch
                    warning('obj.Lick_SideTime() got an error, setting values to NaN');
                    Licking_Side(i_trial,1) = nan;
                    Response_Time(i_trial,1) = nan;
                end
                
            end
            
            obj.licking.time=Licking_Timing;    
            obj.licking.side=Licking_Side;
            obj.licking.response_time = Response_Time;
        end

                
        %% load wavesurfer data (NL added 9/6/12) modified by Guang 12/19/17, by Guang 09/16/2023
        function load_Wavesurfer_Data(obj, wavesurfer_filename)

            
            % establish the sub-fields
            obj.wavesurfer.timestamp = [];
            obj.wavesurfer.trig_trace = [];
            obj.wavesurfer.bitcode_trace = [];
            obj.wavesurfer.pole_trace = [];
            obj.wavesurfer.cue_trace = [];
            obj.wavesurfer.intan_trig_trace = [];
            obj.wavesurfer.photo_input_trace1 = [];
            obj.wavesurfer.photo_input_trace2 = [];
            obj.wavesurfer.photo_input_trace3 = [];
            obj.wavesurfer.video_trigger_input1 = [];
            obj.wavesurfer.video_trigger_input2 = [];
            obj.wavesurfer.video_trigger_input3 = [];
            obj.wavesurfer.video_trigger_input4 = [];

            % load file & parse data
            i_str = findstr(wavesurfer_filename,'.h5');
            wave_file = hdf5info(wavesurfer_filename(1:i_str-1));
            
            
            inverseChannelScales = 1./hdf5read(wave_file.GroupHierarchy.Groups(1).Groups(1).Datasets(3)); % '/header/Acquisition/AnalogChannelScales'
            combinedScaleFactors = 3.0517578125e-4 * inverseChannelScales; % counts-> volts at AI, 3.0517578125e-4 == 10/2^(16-1)
            
            
            dataset_tmp = hdf5read(wave_file.GroupHierarchy.Groups(2).Datasets(1)); % '/sweep_0001/analogScans'
            
            wavesurfer_trig = dataset_tmp(:,1);
            bitcode = dataset_tmp(:,2);
            pole = dataset_tmp(:,3);
            cue = dataset_tmp(:,4);
            intan_trig = dataset_tmp(:,5);
            if size(dataset_tmp,2)>5
                photo_input1 = dataset_tmp(:,6);
            end
            if size(dataset_tmp,2)>6
                photo_input2 = dataset_tmp(:,7);
            end
            if size(dataset_tmp,2)>7
                photo_input3 = dataset_tmp(:,8);
            end
            if size(dataset_tmp,2)>8
            video_trigger_input1 = dataset_tmp(:,9);
            end
            if size(dataset_tmp,2)>9
            video_trigger_input2 = dataset_tmp(:,10);
            end
            if size(dataset_tmp,2)>10
            video_trigger_input3 = dataset_tmp(:,11);
            end
            if size(dataset_tmp,2)>11
            video_trigger_input4 = dataset_tmp(:,12);
            end
            clear dataset_tmp
            
            TimeStamp = (1:size(wavesurfer_trig,1))/obj.wavesurfer_fs;
            
            t_trig_off = [];
            t_pole_on = [];
            t_pole_off = [];
            t_cue_on = [];
            t_cue_off = [];
            

            % segment trials based on triger pulses
            wavesurfer_trig(1)=0;
            wavesurfer_trig(find(wavesurfer_trig<0))=0;
%             i_trig_on = find(diff(wavesurfer_trig)>4000)+1;  % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
%             i_trig_off= find(diff(wavesurfer_trig)<-4000)+1;  
            i_trig_on = find(diff(wavesurfer_trig>6000)==1)+1;
            i_trig_off = find(diff(wavesurfer_trig>6000)==-1)+1;
            idx=find(i_trig_on(2:end)-i_trig_off(1:end-1)>40000);
            i_trig_on=[i_trig_on(1);i_trig_on(idx+1)];
            i_trig_off=[i_trig_off(idx);i_trig_off(end)];
            idx=find(intan_trig(i_trig_on+100)>4000);
            i_trig_on=i_trig_on(idx);
            i_trig_off=i_trig_off(idx);

            disp(['........ obtained ', num2str(min([length(i_trig_on) length(i_trig_off)])), ' trials from Wavesurfer ..........']);
            disp(['processing ... ']);
            
            
            for i_rep = [1:min([length(i_trig_on) length(i_trig_off)])]

%                 if rem(i_rep,10)==0
                    disp([num2str(i_rep), ' trials']);
%                 end
                
                % find segment based on trigger
                t_trig_on_iRep = TimeStamp(i_trig_on(i_rep));
                
                if i_trig_off(1)-i_trig_on(1)<5*obj.wavesurfer_fs
                i_iRep = find(TimeStamp>=t_trig_on_iRep & TimeStamp<(t_trig_on_iRep+6.5));
                if length(i_iRep)>130000
                    i_iRep = i_iRep(1:130000);
                end
                else
                    if length(find(TimeStamp>=t_trig_on_iRep))>=300000
                i_iRep = find(TimeStamp>=t_trig_on_iRep & TimeStamp<(t_trig_on_iRep+15));
                if length(i_iRep)>300000
                    i_iRep = i_iRep(1:300000);
                end
                    else
                        i_iRep = find(TimeStamp>=t_trig_on_iRep);
                    end
                end
              
                
                if size(obj.wavesurfer.timestamp,1)>1 & length(i_iRep)>size(obj.wavesurfer.timestamp,2)
                    warning(['sample size for trial#',num2str(i_rep),' exceeded by ', num2str(length(i_iRep)-size(obj.wavesurfer.timestamp,2)),'; extra samples discarded.']);
                    i_iRep = i_iRep(1:size(obj.wavesurfer.timestamp,2));
                end
                
                TimeStamp_iRep = TimeStamp(i_iRep);                
                TimeStamp_iRep = TimeStamp_iRep-TimeStamp_iRep(1); % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
                
                % get data
                wavesurfer_trig_iRep = double(wavesurfer_trig(i_iRep))*combinedScaleFactors(1);
                bitcode_iRep = double(bitcode(i_iRep))*combinedScaleFactors(2);
                pole_iRep = double(pole(i_iRep))*combinedScaleFactors(3);
                cue_iRep = double(cue(i_iRep))*combinedScaleFactors(4);
                intan_trig_iRep = double(intan_trig(i_iRep))*combinedScaleFactors(5);
                if length(combinedScaleFactors)>5
                    photo_input1_iRep = double(photo_input1(i_iRep))*combinedScaleFactors(6);
                end
                if length(combinedScaleFactors)>6
                    photo_input2_iRep = double(photo_input2(i_iRep))*combinedScaleFactors(7);
                end
                if length(combinedScaleFactors)>7
                    photo_input3_iRep = double(photo_input3(i_iRep))*combinedScaleFactors(8);
                end
                if length(combinedScaleFactors)>8
                    video_trigger_input1_iRep = double(video_trigger_input1(i_iRep))*combinedScaleFactors(9);
                end
                if length(combinedScaleFactors)>9
                    video_trigger_input2_iRep = double(video_trigger_input2(i_iRep))*combinedScaleFactors(10);
                end
                if length(combinedScaleFactors)>10
                    video_trigger_input3_iRep = double(video_trigger_input3(i_iRep))*combinedScaleFactors(11);
                end
                if length(combinedScaleFactors)>11
                    video_trigger_input4_iRep = double(video_trigger_input4(i_iRep))*combinedScaleFactors(12);
                end
                
                if i_trig_off(1)-i_trig_on(1)>5*obj.wavesurfer_fs
                if i_trig_off(i_rep)-i_trig_on(i_rep)<5*obj.wavesurfer_fs
                    if length(TimeStamp_iRep)<300000
                    TimeStamp_iRep(length(TimeStamp_iRep)+1:300000)=[length(TimeStamp_iRep):300000-1]/200000;
                    end
                    wavesurfer_trig_iRep(130000:end)=0;
                    bitcode_iRep(130000:end)=0;
                    pole_iRep(130000:end)=0;
                    cue_iRep(130000:end)=0;
                    intan_trig_iRep(130000:end)=0;
                    if length(combinedScaleFactors)>5
                    photo_input1_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>6
                    photo_input2_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>7
                    photo_input3_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>8
                    video_trigger_input1_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>9
                    video_trigger_input2_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>10
                    video_trigger_input3_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>11
                    video_trigger_input4_iRep(130000:end)=0;
                    end
                end
                end
                
                % single acquisition
                if length(unique(obj.soloTrialIndex(:,1)))==1
                    
                    % get trial # from bitcode
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
%                     if trial_num<i_rep
%                         trial_num=i_rep;
%                     end
                    
                else
                    % get trial # from bitcode
                    if max(bitcode_iRep)>2
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
                    else
                        trial_num=obj.soloTrialIndex(i_rep,2);
                    end
                    
                    % get acq num
%                     if i_rep<=254
                    acq_num = obj.soloTrialIndex(i_rep,1);
%                     else
%                     acq_num = obj.soloTrialIndex(i_rep-1,1);
%                     end
                    
                    trial_num = find(obj.soloTrialIndex(:,1)==acq_num & obj.soloTrialIndex(:,2)==trial_num);
                
                end
                
                
                % save wavesurfer data
                
                obj.wavesurfer.timestamp(trial_num,:) = TimeStamp_iRep;
                obj.wavesurfer.trig_trace(trial_num,:) = wavesurfer_trig_iRep;
                obj.wavesurfer.bitcode_trace(trial_num,:) = bitcode_iRep;
                obj.wavesurfer.pole_trace(trial_num,:) = pole_iRep;
                obj.wavesurfer.cue_trace(trial_num,:) = cue_iRep;
                obj.wavesurfer.intan_trig_trace(trial_num,:) = intan_trig_iRep;
                if length(combinedScaleFactors)>5
                    obj.wavesurfer.photo_input_trace1(trial_num,:) = photo_input1_iRep;
                end
                if length(combinedScaleFactors)>6
                    obj.wavesurfer.photo_input_trace2(trial_num,:) = photo_input2_iRep;
                end
                if length(combinedScaleFactors)>7
                    obj.wavesurfer.photo_input_trace3(trial_num,:) = photo_input3_iRep;
                end
                
                if length(combinedScaleFactors)>8
                    obj.wavesurfer.video_trigger_input1(trial_num,:) = video_trigger_input1_iRep;
                end
                if length(combinedScaleFactors)>9
                    obj.wavesurfer.video_trigger_input2(trial_num,:) = video_trigger_input2_iRep;
                end
                if length(combinedScaleFactors)>10
                    obj.wavesurfer.video_trigger_input3(trial_num,:) = video_trigger_input3_iRep;
                end
                if length(combinedScaleFactors)>11
                    obj.wavesurfer.video_trigger_input4(trial_num,:) = video_trigger_input4_iRep;
                end                 
            end
        end
            
    
        %% load wavesurfer data (NL added 9/6/12) modified by Guang 12/19/17, by Guang 12/27/2025
        function load_Wavesurfer_Data_optobehavior(obj, wavesurfer_filename)

            
            % establish the sub-fields
            obj.wavesurfer.timestamp = [];
            obj.wavesurfer.trig_trace = [];
            obj.wavesurfer.bitcode_trace = [];
            obj.wavesurfer.pole_trace = [];
            obj.wavesurfer.cue_trace = [];
            obj.wavesurfer.intan_trig_trace = [];
            obj.wavesurfer.photo_input_trace1 = [];
            obj.wavesurfer.photo_input_trace2 = [];
            obj.wavesurfer.photo_input_trace3 = [];
            obj.wavesurfer.video_trigger_input1 = [];
            obj.wavesurfer.video_trigger_input2 = [];
            obj.wavesurfer.video_trigger_input3 = [];
            obj.wavesurfer.video_trigger_input4 = [];

            % load file & parse data
            i_str = findstr(wavesurfer_filename,'.h5');
            wave_file = hdf5info(wavesurfer_filename(1:i_str-1));
            
            
            inverseChannelScales = 1./hdf5read(wave_file.GroupHierarchy.Groups(1).Groups(1).Datasets(3)); % '/header/Acquisition/AnalogChannelScales'
            combinedScaleFactors = 3.0517578125e-4 * inverseChannelScales; % counts-> volts at AI, 3.0517578125e-4 == 10/2^(16-1)
            
            
            dataset_tmp = hdf5read(wave_file.GroupHierarchy.Groups(2).Datasets(1)); % '/sweep_0001/analogScans'
            
            wavesurfer_trig = dataset_tmp(:,1);
            bitcode = dataset_tmp(:,2);
            pole = dataset_tmp(:,3);
            cue = dataset_tmp(:,4);
            intan_trig = dataset_tmp(:,5);
            if size(dataset_tmp,2)>5
                photo_input1 = dataset_tmp(:,6);
            end
            if size(dataset_tmp,2)>6
                photo_input2 = dataset_tmp(:,7);
            end
            if size(dataset_tmp,2)>7
                photo_input3 = dataset_tmp(:,8);
            end
            if size(dataset_tmp,2)>8
            video_trigger_input1 = dataset_tmp(:,9);
            end
            if size(dataset_tmp,2)>9
            video_trigger_input2 = dataset_tmp(:,10);
            end
            if size(dataset_tmp,2)>10
            video_trigger_input3 = dataset_tmp(:,11);
            end
            if size(dataset_tmp,2)>11
            video_trigger_input4 = dataset_tmp(:,12);
            end
            clear dataset_tmp
            
            TimeStamp = (1:size(wavesurfer_trig,1))/obj.wavesurfer_fs;
            
            t_trig_off = [];
            t_pole_on = [];
            t_pole_off = [];
            t_cue_on = [];
            t_cue_off = [];
            

            % segment trials based on triger pulses
            wavesurfer_trig(1)=0;
            intan_trig(1)=0;
            wavesurfer_trig(end)=0;
            intan_trig(end)=0;
            wavesurfer_trig(find(wavesurfer_trig<0))=0;
            % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
            i_trig_on = find(diff(wavesurfer_trig>6000)==1)+1;
            i_trig_off = find(diff(wavesurfer_trig>6000)==-1)+1;
            idx=find(i_trig_on(2:end)-i_trig_off(1:end-1)>40000);
            i_trig_on=[i_trig_on(1);i_trig_on(idx+1)];
            i_trig_off=[i_trig_off(idx);i_trig_off(end)];
            idx=find(intan_trig(i_trig_on+100)>4000);
            i_trig_on=i_trig_on(idx);
            i_trig_off=i_trig_off(idx);
            intan_trig_on=find(diff(intan_trig>4000)==1)+1;
            intan_trig_off=find(diff(intan_trig>4000)==-1)+1;

            intan_trig_on=[intan_trig_on(1);intan_trig_on(find(abs(intan_trig_off(1:end-1)-intan_trig_on(2:end))>0.5*obj.wavesurfer_fs)+1)];
            intan_trig_off=[intan_trig_off(find(abs(intan_trig_off(1:end-1)-intan_trig_on(2:end))>0.5*obj.wavesurfer_fs));intan_trig_off(end)];

            trial_length=max(intan_trig_off-intan_trig_on)+10;
            length(i_trig_on)
            length(intan_trig_on)
            length(i_trig_off)
            length(intan_trig_off)

            if length(i_trig_on)~=length(intan_trig_on)|length(i_trig_off)~=length(intan_trig_off)|length(i_trig_on)~=length(i_trig_off)
               warning(['trig on off is not equal']);
            end
            disp(['........ obtained ', num2str(min([length(i_trig_on) length(i_trig_off)])), ' trials from Wavesurfer ..........']);
            disp(['processing ... ']);
            
            
            for i_rep = [1:min([length(i_trig_on) length(i_trig_off)])]

%                 if rem(i_rep,10)==0
                    disp([num2str(i_rep), ' trials']);
%                 end
                
                % find segment based on trigger
                t_trig_on_iRep = TimeStamp(i_trig_on(i_rep));
                t_trig_off_iRep = TimeStamp(intan_trig_off(i_rep));
                
                i_iRep = find(TimeStamp>=t_trig_on_iRep & TimeStamp<=t_trig_off_iRep);
                
                if size(obj.wavesurfer.timestamp,1)>1 & length(i_iRep)>size(obj.wavesurfer.timestamp,2)
                    warning(['sample size for trial#',num2str(i_rep),' exceeded by ', num2str(length(i_iRep)-size(obj.wavesurfer.timestamp,2)),'; extra samples discarded.']);
                    i_iRep = i_iRep(1:size(obj.wavesurfer.timestamp,2));
                end
                
                TimeStamp_iRep = TimeStamp(1:trial_length);                
                TimeStamp_iRep = TimeStamp_iRep-TimeStamp_iRep(1); % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
                
                % get data
                wavesurfer_trig_iRep=zeros(1,trial_length);
                wavesurfer_trig_iRep(1:length(i_iRep)) = double(wavesurfer_trig(i_iRep))*combinedScaleFactors(1);
                bitcode_iRep=zeros(1,trial_length);
                bitcode_iRep(1:length(i_iRep)) = double(bitcode(i_iRep))*combinedScaleFactors(2);
                pole_iRep=zeros(1,trial_length);
                pole_iRep(1:length(i_iRep)) = double(pole(i_iRep))*combinedScaleFactors(3);
                cue_iRep=zeros(1,trial_length);
                cue_iRep(1:length(i_iRep)) = double(cue(i_iRep))*combinedScaleFactors(4);
                intan_trig_iRep=zeros(1,trial_length);
                intan_trig_iRep(1:length(i_iRep)) = double(intan_trig(i_iRep))*combinedScaleFactors(5);
                if length(combinedScaleFactors)>5
                    photo_input1_iRep=zeros(1,trial_length);
                    photo_input1_iRep(1:length(i_iRep)) = double(photo_input1(i_iRep))*combinedScaleFactors(6);
                end
                if length(combinedScaleFactors)>6
                    photo_input2_iRep=zeros(1,trial_length);
                    photo_input2_iRep(1:length(i_iRep)) = double(photo_input2(i_iRep))*combinedScaleFactors(7);
                end
                if length(combinedScaleFactors)>7
                    photo_input3_iRep=zeros(1,trial_length);
                    photo_input3_iRep(1:length(i_iRep)) = double(photo_input3(i_iRep))*combinedScaleFactors(8);
                end
                if length(combinedScaleFactors)>8
                    video_trigger_input1_iRep=zeros(1,trial_length);
                    video_trigger_input1_iRep(1:length(i_iRep)) = double(video_trigger_input1(i_iRep))*combinedScaleFactors(9);
                end
                if length(combinedScaleFactors)>9
                    video_trigger_input2_iRep=zeros(1,trial_length);
                    video_trigger_input2_iRep(1:length(i_iRep)) = double(video_trigger_input2(i_iRep))*combinedScaleFactors(10);
                end
                if length(combinedScaleFactors)>10
                    video_trigger_input3_iRep=zeros(1,trial_length);
                    video_trigger_input3_iRep(1:length(i_iRep)) = double(video_trigger_input3(i_iRep))*combinedScaleFactors(11);
                end
                if length(combinedScaleFactors)>11
                    video_trigger_input4_iRep=zeros(1,trial_length);
                    video_trigger_input4_iRep(1:length(i_iRep)) = double(video_trigger_input4(i_iRep))*combinedScaleFactors(12);
                end
                
                % single acquisition
                if length(unique(obj.soloTrialIndex(:,1)))==1
                    
                    % get trial # from bitcode
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
%                     if trial_num<i_rep
%                         trial_num=i_rep;
%                     end
                    
                else
                    % get trial # from bitcode
                    if max(bitcode_iRep)>2
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
                    else
                        trial_num=obj.soloTrialIndex(i_rep,2);
                    end
                    
                    % get acq num
%                     if i_rep<=254
                    acq_num = obj.soloTrialIndex(i_rep,1);
%                     else
%                     acq_num = obj.soloTrialIndex(i_rep-1,1);
%                     end
                    
                    trial_num = find(obj.soloTrialIndex(:,1)==acq_num & obj.soloTrialIndex(:,2)==trial_num);
                
                end
                
                
                % save wavesurfer data
                
                obj.wavesurfer.timestamp(trial_num,:) = TimeStamp_iRep;
                obj.wavesurfer.trig_trace(trial_num,:) = wavesurfer_trig_iRep;
                obj.wavesurfer.bitcode_trace(trial_num,:) = bitcode_iRep;
                obj.wavesurfer.pole_trace(trial_num,:) = pole_iRep;
                obj.wavesurfer.cue_trace(trial_num,:) = cue_iRep;
                obj.wavesurfer.intan_trig_trace(trial_num,:) = intan_trig_iRep;
                if length(combinedScaleFactors)>5
                    obj.wavesurfer.photo_input_trace1(trial_num,:) = photo_input1_iRep;
                end
                if length(combinedScaleFactors)>6
                    obj.wavesurfer.photo_input_trace2(trial_num,:) = photo_input2_iRep;
                end
                if length(combinedScaleFactors)>7
                    obj.wavesurfer.photo_input_trace3(trial_num,:) = photo_input3_iRep;
                end
                
                if length(combinedScaleFactors)>8
                    obj.wavesurfer.video_trigger_input1(trial_num,:) = video_trigger_input1_iRep;
                end
                if length(combinedScaleFactors)>9
                    obj.wavesurfer.video_trigger_input2(trial_num,:) = video_trigger_input2_iRep;
                end
                if length(combinedScaleFactors)>10
                    obj.wavesurfer.video_trigger_input3(trial_num,:) = video_trigger_input3_iRep;
                end
                if length(combinedScaleFactors)>11
                    obj.wavesurfer.video_trigger_input4(trial_num,:) = video_trigger_input4_iRep;
                end                 
            end
        end
    
    
        %% load wavesurfer data (NL added 5/24/17) modified by Guang 12/19/17, by Guang 03/04/2025
        function load_MultipleWavesurfer_Data(obj, wavesurfer_filename0)

            
            % establish the sub-fields
            obj.wavesurfer.timestamp = [];
            obj.wavesurfer.trig_trace = [];
            obj.wavesurfer.bitcode_trace = [];
            obj.wavesurfer.pole_trace = [];
            obj.wavesurfer.cue_trace = [];
            obj.wavesurfer.intan_trig_trace = [];
            obj.wavesurfer.photo_input_trace1 = [];
            obj.wavesurfer.photo_input_trace2 = [];
            obj.wavesurfer.photo_input_trace3 = [];
            obj.wavesurfer.video_trigger_input1 = [];
            obj.wavesurfer.video_trigger_input2 = [];
            obj.wavesurfer.video_trigger_input3 = [];
            obj.wavesurfer.video_trigger_input4 = [];
            
            for jw=1:length(wavesurfer_filename0)
                wavesurfer_filename=wavesurfer_filename0{jw,1};
            % load file & parse data
            i_str = findstr(wavesurfer_filename,'.h5');
            wave_file = hdf5info(wavesurfer_filename(1:i_str-1));
            
            
            inverseChannelScales = 1./hdf5read(wave_file.GroupHierarchy.Groups(1).Groups(1).Datasets(3)); % '/header/Acquisition/AnalogChannelScales'
            combinedScaleFactors = 3.0517578125e-4 * inverseChannelScales; % counts-> volts at AI, 3.0517578125e-4 == 10/2^(16-1)
            
            
            dataset_tmp = hdf5read(wave_file.GroupHierarchy.Groups(2).Datasets(1)); % '/sweep_0001/analogScans'
            
            wavesurfer_trig = dataset_tmp(:,1);
            bitcode = dataset_tmp(:,2);
            pole = dataset_tmp(:,3);
            cue = dataset_tmp(:,4);
            intan_trig = dataset_tmp(:,5);
            if size(dataset_tmp,2)>5
                photo_input1 = dataset_tmp(:,6);
            end
            if size(dataset_tmp,2)>6
                photo_input2 = dataset_tmp(:,7);
            end
            if size(dataset_tmp,2)>7
                photo_input3 = dataset_tmp(:,8);
            end          
            if size(dataset_tmp,2)>8
            video_trigger_input1 = dataset_tmp(:,9);
            end
            if size(dataset_tmp,2)>9
            video_trigger_input2 = dataset_tmp(:,10);
            end
            if size(dataset_tmp,2)>10
            video_trigger_input3 = dataset_tmp(:,11);
            end
            if size(dataset_tmp,2)>11
            video_trigger_input4 = dataset_tmp(:,12);
            end
            clear dataset_tmp
            
            TimeStamp = (1:size(wavesurfer_trig,1))/obj.wavesurfer_fs;
            
            t_trig_off = [];
            t_pole_on = [];
            t_pole_off = [];
            t_cue_on = [];
            t_cue_off = [];
            

            % segment trials based on triger pulses
            wavesurfer_trig(1)=0;
            wavesurfer_trig(find(wavesurfer_trig<0))=0;
%             i_trig_on = find(diff(wavesurfer_trig)>4000)+1;  % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
%             i_trig_off= find(diff(wavesurfer_trig)<-4000)+1;  
            i_trig_on = find(diff(wavesurfer_trig>6000)==1)+1;
            i_trig_off = find(diff(wavesurfer_trig>6000)==-1)+1;
            idx=find(i_trig_on(2:end)-i_trig_off(1:end-1)>40000);
            i_trig_on=[i_trig_on(1);i_trig_on(idx+1)];
            i_trig_off=[i_trig_off(idx);i_trig_off(end)];
            idx=find(intan_trig(i_trig_on+100)>4000);
            i_trig_on=i_trig_on(idx);
            i_trig_off=i_trig_off(idx);

            disp(['........ obtained ', num2str(min([length(i_trig_on) length(i_trig_off)])), ' trials from Wavesurfer ..........']);
            disp(['processing ... ']);
            
            
            for i_rep = [1:min([length(i_trig_on) length(i_trig_off)])]

%                 if rem(i_rep,10)==0
                    disp([num2str(i_rep), ' trials']);
%                 end
                
                % find segment based on trigger
                t_trig_on_iRep = TimeStamp(i_trig_on(i_rep));
                
                if i_trig_off(1)-i_trig_on(1)<5*obj.wavesurfer_fs
                i_iRep = find(TimeStamp>=t_trig_on_iRep & TimeStamp<(t_trig_on_iRep+6.5));
                if length(i_iRep)>130000
                    i_iRep = i_iRep(1:130000);
                end
                else
                    if length(find(TimeStamp>=t_trig_on_iRep))>=300000
                i_iRep = find(TimeStamp>=t_trig_on_iRep & TimeStamp<(t_trig_on_iRep+15));
                if length(i_iRep)>300000
                    i_iRep = i_iRep(1:300000);
                end
                    else
                        i_iRep = find(TimeStamp>=t_trig_on_iRep);
                    end
                end
              
                
                if size(obj.wavesurfer.timestamp,1)>1 & length(i_iRep)>size(obj.wavesurfer.timestamp,2)
                    warning(['sample size for trial#',num2str(i_rep),' exceeded by ', num2str(length(i_iRep)-size(obj.wavesurfer.timestamp,2)),'; extra samples discarded.']);
                    i_iRep = i_iRep(1:size(obj.wavesurfer.timestamp,2));
                end
                
                TimeStamp_iRep = TimeStamp(i_iRep);                
                TimeStamp_iRep = TimeStamp_iRep-TimeStamp_iRep(1); % the first sample is the first sample after the trigger ON, so at t=0, trigger is already ON.
                
                % get data
                wavesurfer_trig_iRep = double(wavesurfer_trig(i_iRep))*combinedScaleFactors(1);
                bitcode_iRep = double(bitcode(i_iRep))*combinedScaleFactors(2);
                pole_iRep = double(pole(i_iRep))*combinedScaleFactors(3);
                cue_iRep = double(cue(i_iRep))*combinedScaleFactors(4);
                intan_trig_iRep = double(intan_trig(i_iRep))*combinedScaleFactors(5);
                if length(combinedScaleFactors)>5
                    photo_input1_iRep = double(photo_input1(i_iRep))*combinedScaleFactors(6);
                end
                if length(combinedScaleFactors)>6
                    photo_input2_iRep = double(photo_input2(i_iRep))*combinedScaleFactors(7);
                end
                if length(combinedScaleFactors)>7
                    photo_input3_iRep = double(photo_input3(i_iRep))*combinedScaleFactors(8);
                end
                if length(combinedScaleFactors)>8
                    video_trigger_input1_iRep = double(video_trigger_input1(i_iRep))*combinedScaleFactors(9);
                end
                if length(combinedScaleFactors)>9
                    video_trigger_input2_iRep = double(video_trigger_input2(i_iRep))*combinedScaleFactors(10);
                end
                if length(combinedScaleFactors)>10
                    video_trigger_input3_iRep = double(video_trigger_input3(i_iRep))*combinedScaleFactors(11);
                end
                if length(combinedScaleFactors)>11
                    video_trigger_input4_iRep = double(video_trigger_input4(i_iRep))*combinedScaleFactors(12);
                end
                
                if i_trig_off(1)-i_trig_on(1)>5*obj.wavesurfer_fs
                if i_trig_off(i_rep)-i_trig_on(i_rep)<5*obj.wavesurfer_fs
                    if length(TimeStamp_iRep)<300000
                    TimeStamp_iRep(length(TimeStamp_iRep)+1:300000)=[length(TimeStamp_iRep):300000-1]/200000;
                    end
                    wavesurfer_trig_iRep(130000:end)=0;
                    bitcode_iRep(130000:end)=0;
                    pole_iRep(130000:end)=0;
                    cue_iRep(130000:end)=0;
                    intan_trig_iRep(130000:end)=0;
                    if length(combinedScaleFactors)>5
                    photo_input1_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>6
                    photo_input2_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>7
                    photo_input3_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>8
                    video_trigger_input1_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>9
                    video_trigger_input2_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>10
                    video_trigger_input3_iRep(130000:end)=0;
                    end
                    if length(combinedScaleFactors)>11
                    video_trigger_input4_iRep(130000:end)=0;
                    end
                end
                end
                
                % single acquisition
                if length(unique(obj.soloTrialIndex(:,1)))==1
                    
                    % get trial # from bitcode
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
%                     if trial_num<i_rep
%                         trial_num=i_rep;
%                     end
                    
                else
                    % get trial # from bitcode
                    if max(bitcode_iRep)>2
                    trial_num = obj.func_read_bitcode(bitcode_iRep,TimeStamp_iRep);
                    else
                        trial_num=obj.soloTrialIndex(i_rep,2);
                    end
                    % get acq num
%                     if i_rep<=254
                    acq_num = obj.soloTrialIndex(i_rep,1);
%                     else
%                     acq_num = obj.soloTrialIndex(i_rep-1,1);
%                     end
                    
                    trial_num = find(obj.soloTrialIndex(:,1)==acq_num & obj.soloTrialIndex(:,2)==trial_num);
                
                end
                
                
                % save wavesurfer data
                
                obj.wavesurfer.timestamp(trial_num,:) = TimeStamp_iRep;
                obj.wavesurfer.trig_trace(trial_num,:) = wavesurfer_trig_iRep;
                obj.wavesurfer.bitcode_trace(trial_num,:) = bitcode_iRep;
                obj.wavesurfer.pole_trace(trial_num,:) = pole_iRep;
                obj.wavesurfer.cue_trace(trial_num,:) = cue_iRep;
                obj.wavesurfer.intan_trig_trace(trial_num,:) = intan_trig_iRep;
                if length(combinedScaleFactors)>5
                    obj.wavesurfer.photo_input_trace1(trial_num,:) = photo_input1_iRep;
                end
                if length(combinedScaleFactors)>6
                    obj.wavesurfer.photo_input_trace2(trial_num,:) = photo_input2_iRep;
                end
                if length(combinedScaleFactors)>7
                    obj.wavesurfer.photo_input_trace3(trial_num,:) = photo_input3_iRep;
                end
                
                if length(combinedScaleFactors)>8
                    obj.wavesurfer.video_trigger_input1(trial_num,:) = video_trigger_input1_iRep;
                end
                if length(combinedScaleFactors)>9
                    obj.wavesurfer.video_trigger_input2(trial_num,:) = video_trigger_input2_iRep;
                end
                if length(combinedScaleFactors)>10
                    obj.wavesurfer.video_trigger_input3(trial_num,:) = video_trigger_input3_iRep;
                end
                if length(combinedScaleFactors)>11
                    obj.wavesurfer.video_trigger_input4(trial_num,:) = video_trigger_input4_iRep;
                end               
            end
            end
        end
            
                
            
        
        
        %% a function to read out trial number from bitcode (NL added 9/6/12)        
        function trial_NO = func_read_bitcode(obj, bitcode, time_stamp)
            
            
            time_stamp = time_stamp*1000;  % input should be in sec, convert to ms here
            

            threshold = 0.8;
            bitcode(1) = 0;
            bitcode(end) = 0;
            i_bitcode_on = find(diff(bitcode>max(double(bitcode))*threshold)==1);
            i_bitcode_off = find(diff(bitcode>max(double(bitcode))*threshold)==-1);
            
            idx=find(bitcode>max(double(bitcode))*threshold);

            if time_stamp(idx(1))>5000
            i_start = 10000.5;
            elseif time_stamp(idx(1))>2000
                    i_start = 3500.5;
            elseif time_stamp(idx(1))>1000
                    i_start = 1500.5;
            else
                i_start=500.5;
            end
            
            bitcode_Interval = 7;  % samples       % 2 ms for bit, 5 ms for gap
            numBit = 10;
            t_start = (0:numBit-1)*bitcode_Interval + i_start;
            t_end = t_start+1.5;
%             t_end=[t_start(1:end-1)+1.5 t_start(end)+50];%trial no. 512 bit pulse occassionally happened later
            bit_tmp = [];
            

            for i_bit = 1:numBit
                
                i_sample = find(time_stamp>=t_start(i_bit) & time_stamp<=t_end(i_bit));
%                 bit_state = mean(bitcode(i_sample));

                bit_state = max(bitcode(i_sample));%trial no. 512 bit pulse occassionally happened later
                
                if bit_state>max(double(bitcode))*threshold
                    bit_tmp = [bit_tmp '1'];
                else
                    bit_tmp = [bit_tmp '0'];
                end

            end
            i_str = findstr(bit_tmp,'1');
            bit_tmp = fliplr(bit_tmp(1:i_str(end)));
            
            trial_NO = bin2dec(bit_tmp);


        end
            

        
             
%         function Unit( obj, single_unit_dir )
%             
%             %%%%%%%%%%%%%%%%%%% NL 9/23/12 %%%%%%%%%%%%%
%             %% original
%             %             single_unit_filelist = dir([single_unit_dir,'SingleUnit*']);
%             %
%             %             % go through the sorted units one by one
%             %             for i_file = 1:size(single_unit_filelist,1)
%             %
%             %                 disp(['--------------------------------------------']);
%             %                 disp(['processing file ',single_unit_filelist(i_file).name]);
%             %
%             %                 % load the data
%             %                 load([single_unit_dir, single_unit_filelist(i_file).name]);
%             %                 obj.units{i_file}=unit;
%             %
%             %             end
%             
% 
%             % the same as old code
%             single_unit_filelist = dir([single_unit_dir,'SingleUnit*']);
%             
%             % go through the sorted units one by one
%             for i_file = 1:size(single_unit_filelist,1)
%                 
%                 disp(['--------------------------------------------']);
%                 disp(['processing file ',single_unit_filelist(i_file).name]);
%                 
%                 % load the data
%                 load([single_unit_dir, single_unit_filelist(i_file).name]);
%                 obj.units{i_file}=unit;
%                 
%             end
%             
%             
%         end
        
        
        
    end
    methods( Hidden = true )
        
%         function newSaved=Change_Saved_Field_Names(obj, saved)
%             
%             fieldNames = fieldnames(saved);
%             protocolName=[];
%             for i=1:size(fieldNames, 1)
%                 if( strfind(fieldNames{i}, '_hit_history') )
%                     index=strfind(fieldNames{i}, '_hit_history');
%                     protocolName=fieldNames{i}(1:index-1);
%                 end
%             end
%             
%             if( isempty(protocolName) )
%                 error('didn\''t find protocol name')
%             else
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   add new fields
%                 newSaved=saved;
%                 for i=1:size(fieldNames, 1)
%                     if( strfind(fieldNames{i}, protocolName) )
%                         index=strfind(fieldNames{i}, protocolName);
%                         eval(['newSaved.' fieldNames{i}(1:index-1) 'stim_timing2AFCobj' fieldNames{i}(index+length(protocolName):end) '= saved.' fieldNames{i},';']);
%                     end
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  remove old fields
%                 fieldNames = fieldnames(newSaved);
%                 index_to_remove=[];
%                 for i=1:size(fieldNames, 1)
%                     
%                     if( strfind(fieldNames{i}, protocolName) )
%                         index_to_remove=[index_to_remove i];
%                     end
%                     
%                 end
%                 newSaved=rmfield(newSaved, {fieldNames{index_to_remove}});
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  remove old fields
%             end
%         end
        
         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%<--------this needs to be modified and added for stim protocols
%         function Stim_State_History(obj)
%             
%             sidesTypes=obj.sidesTypes;
%             
%             for i=1:length(sidesTypes)
%                 
%                 if ~isempty(strfind(sidesTypes{i},'sample_period'))
%                     obj.stimStateHistory(i) = 3;
%                 elseif ~isempty(strfind(sidesTypes{i},'delay_period'))
%                     obj.stimStateHistory(i) = 4;
%                 elseif ~isempty(strfind(sidesTypes{i},'reward_period')) | ~isempty(strfind(sidesTypes{i},'response_period'))  %% NL 4/3/13
%                     obj.stimStateHistory(i) = 5;
%                 elseif ~isempty(strfind(sidesTypes{i},'on_whisker_contact'))
%                     obj.stimStateHistory(i) = 0;
%                 %elseif ~isempty(strfind(sidesTypes{i},'r_s')) & isempty(strfind(sidesTypes{i},'on_whisker_contact')) & isempty(strfind(sidesTypes{i},'_period'))  %% NL 4/17/13
%                 %    obj.stimStateHistory(i) = 5;
%                 elseif ~isempty(strfind(sidesTypes{i},'antidromic_tagging_PT')) %% NL 8/5/13
%                     obj.stimStateHistory(i) = 6;
%                 elseif ~isempty(strfind(sidesTypes{i},'antidromic_tagging_IT')) %% NL 8/5/13
%                     obj.stimStateHistory(i) = 7;
%                 elseif ~isempty(strfind(sidesTypes{i},'short_stim_s_late')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 8;
%                 elseif ~isempty(strfind(sidesTypes{i},'short_stim_d_early')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 9;
%                 elseif ~isempty(strfind(sidesTypes{i},'short_stim_d_late')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 10;
%                 elseif ~isempty(strfind(sidesTypes{i},'tagging_IT_s_late')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 11;
%                 elseif ~isempty(strfind(sidesTypes{i},'tagging_IT_d_early')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 12;
%                 elseif ~isempty(strfind(sidesTypes{i},'tagging_IT_d_late')) %% NL 10/3/13
%                     obj.stimStateHistory(i) = 13;
%                 elseif ~isempty(strfind(sidesTypes{i},'on_whisker_contact'))
%                     obj.stimStateHistory(i) = 0;
%                 elseif ~isempty(strfind(sidesTypes{i},'_n'))
%                     % this is just a dummy variable, not actually used
%                     obj.stimStateHistory(i) = 1;
%                 elseif ~isempty(strfind(sidesTypes{i},'r_s')) %% NL 8/5/13
%                     obj.stimStateHistory(i) = 2;        % unspecific stim type (from StimPulse protocol)
%                 elseif ~isempty(strfind(sidesTypes{i},'400d'))
%                     obj.stimStateHistory(i) = 15;
%                 else
%                     error('unrecognized stim type selection')
%                 end
%                 
%             end
%             obj.stimStateHistory=obj.stimStateHistory';
%         end
        
        
        
        function Trim_Trials(obj, trimTrials)
            % tris function trim unnecessary trials
            obj.trials.hitHistory(trimTrials)=0;
            obj.trials.missHistory(trimTrials)=0;
            obj.trials.noResponseHistory(trimTrials)=0;
            obj.trials.EarlyLicksHistory(trimTrials)=0;
            
%             obj.aomPowerHistory(trimTrials)=0;
            %%-----need to add---->obj.xGalvoHistory=saved.stim_timing2AFCobj_stim_xGalvo_pos_history;
            %%-----need to add---->obj.yGalvoHistory=saved.stim_timing2AFCobj_stim_yGalvo_pos_history;
            %%-----need to add---->obj.Stim_State_History;

                
        end
        
%         function Eval_Expression=Trial_Nums(obj, which_trials)
%             
%             %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%%
%             %Eval_Expression=['obj.trials' which_trials '.hitHistory(trimTrials)=0;'];
%             %Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.missHistory(trimTrials)=0;'];
%             %Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.noResponseHistory(trimTrials)=0;'];
%             %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%%
%             
%             Eval_Expression=['obj.trials' which_trials '.hitTrialNums = find(obj.trials' which_trials '.hitHistory==1);'];
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.missTrialNums = find(obj.trials' which_trials '.missHistory==1);'];
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.noResponseTrialNums = find(obj.trials' which_trials '.noResponseHistory==1);'];
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.trialNums =union( union(obj.trials' which_trials '.hitTrialNums, obj.trials' which_trials '.missTrialNums),'...
%                 'obj.trials' which_trials '.noResponseTrialNums);'];
%             % length
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.length =length(obj.trials' which_trials '.trialNums);'];
%             % calculate rate within each category
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.hitRate =length(obj.trials' which_trials '.hitTrialNums)/obj.trials' which_trials '.length;'];
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.missRate =length(obj.trials' which_trials '.missTrialNums)/obj.trials' which_trials '.length;'];
%             Eval_Expression=[Eval_Expression 'obj.trials' which_trials '.noResponseRate =length(obj.trials' which_trials '.noResponseTrialNums)/obj.trials' which_trials '.length;'];
%             
%             %             obj.trialsL.hitHistory(trimTrials)=0;
%             %             obj.trialsL.missHistory(trimTrials)=0;
%             %             obj.trialsL.noResponseHistory(trimTrials)=0;
%             %             obj.trialsL.hitTrialNums = find(obj.trialsL.hitHistory==1);
%             %             obj.trialsL.missTrialNums = find(obj.trialsL.missHistory==1);
%             %             obj.trialsL.noResponseTrialNums=find(obj.trialsL.noResponseHistory==1);
%             %             obj.trialsL.trialNums=union( union(obj.trialsL.hitTrialNums, obj.trialsL.missTrialNums), obj.trialsL.noResponseTrialNums);
%         end
        
%         function Alarm_Which_Trials(obj, which_trials)
%             
%             eval( ['hitTrialNums=obj.trials' which_trials '.hitTrialNums;'] );
%             eval( ['missTrialNums=obj.trials' which_trials '.missTrialNums;'] );
%             eval( ['noResponseTrialNums=obj.trials' which_trials '.noResponseTrialNums;'] );
%             eval( ['len=obj.trials' which_trials '.length;'] );     %the length of total trial nums with each category, including hit,miss, no response
%             
%             %%% calculate alarm nums in hit, miss or no response trials
%             [hitAlarmNums ~]=obj.Alarm_Analysis(hitTrialNums);
%             [missAlarmNums ~]=obj.Alarm_Analysis( missTrialNums );
%             [noResponseAlarmNums ~]=obj.Alarm_Analysis( noResponseTrialNums );
%             
%             hitAlarmRate=length(hitAlarmNums)/length(hitTrialNums);
%             missAlarmRate=length(missAlarmNums)/length(missTrialNums);
%             noResponseAlarmRate=length(noResponseAlarmNums)/length(noResponseTrialNums);
%             alarmNums=union( union(hitAlarmNums, missAlarmNums), noResponseAlarmNums);
%             alarmLen=length(alarmNums);
%             alarmRate=alarmLen/len;
%             
%             %%% calculate in which state the alarm was triggered
%             [~, alarmState]=obj.Alarm_Analysis( alarmNums );
%             alarmNumsSample=alarmState.state53;
%             alarmNumsDelay=setdiff(alarmState.state55, union(union(alarmState.state53,alarmState.state59), alarmState.state61) );
%             alarmNumsPoleUp=setdiff(union(alarmState.state59, alarmState.state61), alarmState.state53 );
%             
%             eval( ['obj.trials' which_trials '.hitAlarmNums=hitAlarmNums' ] );
%             eval( ['obj.trials' which_trials '.missAlarmNums=missAlarmNums' ] );
%             eval( ['obj.trials' which_trials '.noResponseAlarmNums=noResponseAlarmNums' ] );
%             eval( ['obj.trials' which_trials '.alarmNums=alarmNums' ] );
%             eval( ['obj.trials' which_trials '.hitAlarmRate=hitAlarmRate' ] );
%             eval( ['obj.trials' which_trials '.missAlarmRate=missAlarmRate' ] );
%             eval( ['obj.trials' which_trials '.noResponseAlarmRate=noResponseAlarmRate' ] );
%             eval( ['obj.trials' which_trials '.alarmLen=alarmLen' ] );
%             eval( ['obj.trials' which_trials '.alarmRate=alarmRate' ] );
%             eval( ['obj.trials' which_trials '.alarmNumsSample=alarmNumsSample' ] );
%             eval( ['obj.trials' which_trials '.alarmNumsDelay=alarmNumsDelay' ] );
%             eval( ['obj.trials' which_trials '.alarmNumsPoleUp=alarmNumsPoleUp' ] );
%             
%         end
        
%         function [AlarmTrials, AlarmStateTrials]=Alarm_Analysis(obj, trialNums )
%             
%             % the function identifies in which trials alarm was triggered by animal's licking during either sample, delay period or pole up period
%             % AlarmTrials: all the alarm trials within TrialNums
%             % AlarmStateTrials: alarm trials within which a specific state is visited
%             % TrialNums: a vector contains trial numbers, i.e. trialsRNoStim.HitTrialNums, LHitTrialNums
%             % StateMatrixTrialEvents: record of state visited
%             % StateAlarm: vector contains states within which alarm is triggered.
%             % 9/23/2011 by Zengcai Guo
%             
%             if( length(obj.stateAlarm) ==0)
%                 disp( 'must present an alarm state' )
%                 return;
%             end
%             
%             if ( length(trialNums) ==0 )
%                 AlarmTrials=[];
%                 for nstate=1:length(obj.stateAlarm)
%                     eval(['AlarmStateTrials.state' num2str(obj.stateAlarm(nstate)) ' = [];']);
%                 end
%                 return;
%             end
%             
%             
%             AlarmTrials=[];     AlarmStateTrials=[];
%             SM_States=obj.eventsHistory;
%             for nstate=1:length(obj.stateAlarm)
%                 nhit=0;
%                 eval( ['AlarmStateTrials.state' num2str(obj.stateAlarm(nstate)) '=[];'] );
%                 for j=1:length(trialNums)
%                     if( length( find( SM_States{ trialNums(j) } (:,1)==obj.stateAlarm(nstate ) )) )
%                         nhit=nhit+1;
%                         eval(['AlarmStateTrials.state' num2str(obj.stateAlarm(nstate)) '(nhit) = trialNums(j);']);
%                     end
%                 end
%                 eval(['AlarmTrials=union(AlarmTrials, AlarmStateTrials.state' num2str(obj.stateAlarm(nstate)) ');'] );
%             end
%         end
%         
        
        
        
        
        
        
        
        %%%%%%%%%%% NL change (v120912) %%%%%%%%%%%%%%% GC change 20250219
        function addYesNoDiscriminationCIBRData(obj, obj_tmp)
  
%             try
            %%%%%%%%%%% NL change (v121113) %%%%%%%%%%%%%%%
            obj.trim             = cat(1, obj.trim, obj_tmp.trim);
            
            obj.mouseName        = cat(1, obj.mouseName, obj_tmp.mouseName);
            obj.sessionName      = cat(2, obj.sessionName, obj_tmp.sessionName);
            obj.sessionType      = cat(1, obj.sessionType, obj_tmp.sessionType);
            obj.sessionNameTag   = cat(1, obj.sessionNameTag, obj_tmp.sessionNameTag);
            obj.eventsHistory    = cat(1, obj.eventsHistory, obj_tmp.eventsHistory);
            obj.stimType         = cat(1, obj.stimType, obj_tmp.stimType);
            obj.stimProb         = cat(1, obj.stimProb, obj_tmp.stimProb);
            obj.sides            = cat(1, obj.sides, obj_tmp.sides);
            obj.sidesTypes       = cat(1, obj.sidesTypes, obj_tmp.sidesTypes);
            obj.PoleSideType       = cat(1, obj.PoleSideType, obj_tmp.PoleSideType);
            obj.contextType       = cat(1, obj.contextType, obj_tmp.contextType);
            obj.TrialPoleRight       = cat(1, obj.TrialPoleRight, obj_tmp.TrialPoleRight);
            obj.TrialPoleLeft       = cat(1, obj.TrialPoleLeft, obj_tmp.TrialPoleLeft);
            obj.tasteType       = cat(1, obj.tasteType, obj_tmp.tasteType);

            obj.photo_input1History  = cat(1, obj.photo_input1History, obj_tmp.photo_input1History);
            obj.photo_input2History    = cat(1, obj.photo_input2History, obj_tmp.photo_input2History);
            obj.photo_input3History    = cat(1, obj.photo_input3History, obj_tmp.photo_input3History);
            obj.stimStateHistory = cat(1, obj.stimStateHistory, obj_tmp.stimStateHistory);
            
            obj.TrialStartTimestamp=cat(1, obj.TrialStartTimestamp, obj_tmp.TrialStartTimestamp);
            obj.TrialEndTimestamp=cat(1, obj.TrialEndTimestamp, obj_tmp.TrialEndTimestamp);
            
            obj.trials.hitHistory = cat(1, obj.trials.hitHistory, obj_tmp.trials.hitHistory);
            obj.trials.missHistory = cat(1, obj.trials.missHistory, obj_tmp.trials.missHistory);
            obj.trials.noResponseHistory = cat(1, obj.trials.noResponseHistory, obj_tmp.trials.noResponseHistory);
            obj.trials.EarlyLicksHistory = cat(1, obj.trials.EarlyLicksHistory, obj_tmp.trials.EarlyLicksHistory);
            obj.trials.trialNums = cat(1, obj.trials.trialNums, obj_tmp.trials.trialNums);
%             catch
%                 keyboard
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% process trial statistics: hit, miss, no response trial nums, rate, length
            

            obj_tmp.soloTrialIndex(:,1) = obj.soloTrialIndex(end,1)+1;
            obj.soloTrialIndex = cat(1, obj.soloTrialIndex, obj_tmp.soloTrialIndex);
            
            
        end
        
        end
    
end




