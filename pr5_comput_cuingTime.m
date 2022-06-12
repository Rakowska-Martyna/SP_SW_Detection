%% Compute the total duration of cue and no-cue time
clc
clear

%% Define path
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueingPeriod';%output directory
ch_hypnPath	= '/yourDirectory/Data/Sleep_scoring'; %hypnograms direcotry
ch_cuesPath = '/yourDirectory/Data/Sleep_EEG'; %markers/EEG directory


%% Define sleep stage
sleep_stage = [2,3];

        if sleep_stage == [2,3]
            stage = 23;         
        elseif sleep_stage == [2]
            stage = 2;
        elseif sleep_stage == [3]
            stage = 3;
        end
 
 %% Select participants
 % Need to analyse data from each system separately       
 downsampled = 1;  
 
 if downsampled == 1
    mx_pFiles	= [{'MRI_part7_sleep2'},{'MRI_part15_sleep1'}];
    
elseif downsampled == 0
    mx_pFiles	= [{'MRI_part20_sleep1'}];
                
 end

nm_ppnt = numel(mx_pFiles);

%% Define cue and no-cue time in samples
cue_period = 875; % = 3.5 seconds * 250 Hz (sampling rate after downsampling)
break_period = 5000; % 20 * 250Hz

for pp = 1:size(mx_pFiles,2)   
 %% Load hypnogram
        fprintf('Analysing ppnt %s: \n',mx_pFiles{pp})
        fprintf('Loading data hypnogram: ')
        tic
        ch_hypfile  = sprintf('psgHypno-%s.mat',mx_pFiles{pp});
        st_hyp      = load(fullfile(ch_hypnPath,ch_hypfile));
        vt_dHypno   = single(st_hyp.dat(1,:));
        vt_hArous   = st_hyp.arousals{1,1};
        toc
        
        vt_hArous_samples = vt_hArous * 250; % convert seconds to samples      
        arousal.start = vt_hArous_samples(:,1); 
        arousal.end = vt_hArous_samples(:,2);
        
        %% Load cues
        fprintf('Loading cue markers: ')
        tic
        ch_cuefile     = sprintf('%s.vmrk',mx_pFiles{pp});       
        ch_cuefilename = fullfile(ch_cuesPath,ch_cuefile);
        markers = importdata(ch_cuefilename,'',10000);
        cues = markers(11:end,:);  % marker 1 starts at line 11 for this data
        % Extract sample during which the cue was played
        for k = 1:numel(cues)
            x = cues(k,1);
            y = strsplit(x{1,1},',');
            t = str2double(y(:,3));
            vt_cues(k,1) = t;
        end
       
        if downsampled == 1
        vt_cues = vt_cues/2; 
        end              
        
        
 %% Prepare variables       
 stageID = st_hyp.dat';
 timeID = st_hyp.timeEpoch';
 timeIDsamples = timeID * 250;
 correct_epochs = ismember(stageID,sleep_stage); % identify epochs of sleep_stage of interest
 correct_times = timeIDsamples(correct_epochs); % identify 'start' samples of each epoch of interest
 
 nm_cues    = numel(vt_cues); % number of cues
 nm_epochs = numel(correct_times); % number of epochs
 nm_ar = length(vt_hArous_samples); % number of arousals
 
 sto_p = zeros(nm_cues,1);%01 vector where 1 = individual cue (not immediatelly followed by another)
 period = zeros(nm_cues,1);%duration of the cue period   
 period_nocue = zeros(nm_cues,1);%duration of the no-cue period
 
 cue_tosubtract(1,1)  = 0; % if it needs to be removed from cue period
 cue_tosubtract(1,2)  = 0; % how many seconds need to be removed from cue period
 cue_tosubtract(1,3)  = 0; % which cue it is
 cue_tosubtract(1,4)  = 0; % which arousal it is
 
 nocue_tosubtract(1,1)  = 0; % if it needs to be removed from nocue period
 nocue_tosubtract(1,2)  = 0; % how many seconds need to be removed from nocue period
 nocue_tosubtract(1,3)  = 0; % which cue it is
 nocue_tosubtract(1,4)  = 0; % which arousal it is
 
 %% Identify cue and no-cue period
 
  for a = 1:nm_cues % for all cues.. / check each cue...
            if a~=nm_cues % ... but the last one /... except the last one
                for f = 1:nm_epochs % for all correct epochs (i.e. epochs in the correct sleep stage)/ if any of the correct epochs
                    epochduration = correct_times(f,1) + (30 * 250); % next epoch start
                    x = correct_times(f,1):epochduration; % all the samples from the current epoch to the next / (any sample within that epoch)
                     for j = 1:numel(x) % for all the samples from the current epoch to the next
                         % get cue_period = any period between cue and cue
                         % +3.5 s (or if cue is not followed by another one, 
                         % count it as a single cue and add to the total number)
                        
                        if vt_cues(a,1) + cue_period >= vt_cues(a+1,1) && vt_cues(a,1) >= x(j) && vt_cues(a,1) <= x(j)+1 %&& vt_cues(a,1) == x(j)  %if a cue + cue_period is less than the next cue, then the cue period is the time between the 2 cues
                                period(a,1) = vt_cues(a+1,1) - vt_cues(a,1); % duration of cue period
                        elseif vt_cues(a,1) >= x(j) && vt_cues(a,1) <= x(j)+1 %vt_cues(a,1) == x(j) % but if the cue is not immediatelly followed by the next cue (within 3.5s) then the cue period is 3.5s
                                sto_p(a,1) = 1; % it's a cue period        
                        end 
                        
                        % get nocue_period
                        if (vt_cues(a+1,1) - vt_cues(a,1)) >= break_period ...
                        && (vt_cues(a+1,1) - vt_cues(a,1)) <= (break_period + cue_period)...
                        && vt_cues(a,1) >= x(j) && vt_cues(a,1) <= x(j)+1 %&& vt_cues(a+1,1) == x(j) % only the first cue (before the break) within the correct period
                        period_nocue(a,1) = vt_cues(a+1,1) - (vt_cues(a,1)+cue_period);
                        end
                     end
                end
                
%% Remove arousals 

        % remove arousals from cue period
                        if sto_p(a,1) > 0 % if it's an individual cue
                            cuesp = vt_cues(a,1):(vt_cues(a,1)+cue_period);%all samples between 2 cues
                            for h = 1:length(vt_hArous_samples)
                                ar_dur = arousal.start(h,1):arousal.end(h,1);%all samples within an arousal
                                cue_ar = sum(ismember(ar_dur,cuesp));%if a between-cues sample is also an arousal sample then this sample needs to be removed
                                if cue_ar >=1
                                           cue_tosubtract(end+1,1) = 1;
                                           cue_tosubtract(end,2) = cue_ar/250; %in seconds
                                           cue_tosubtract(end,3) = a;
                                           cue_tosubtract(end,4) = h;
                                end 
                            end
                        end
                        
                       if period(a,1) > 0 % if its a cue followed by another cue within 3.5s
                            cuesp_p = vt_cues(a,1):(vt_cues(a,1)+period(a,1));%all samples between 2 cues
                            for h = 1:length(vt_hArous_samples)
                                ar_dur_p = arousal.start(h,1):arousal.end(h,1);%all samples within an arousal
                                cue_ar_p = sum(ismember(ar_dur_p,cuesp_p));%if a between-cues sample is also an arousal sample then this sample needs to be removed
                                if cue_ar_p >=1
                                           cue_tosubtract(end+1,1) = 1;
                                           cue_tosubtract(end,2) = cue_ar_p/250; %in seconds
                                           cue_tosubtract(end,3) = a;
                                           cue_tosubtract(end,4) = h;
                                end 
                            end
                       end
                        
                    % Remove arousals from no-cue period    
                     if period_nocue(a,1) >0 
                            nocuesp = vt_cues(a,1):(vt_cues(a,1)+period_nocue(a,1));
                            for h = 1:length(vt_hArous_samples)
                                ar_dur_no = arousal.start(h,1):arousal.end(h,1);
                                cue_ar_no = sum(ismember(ar_dur_no,nocuesp));
                                if cue_ar_no >=1
                                           nocue_tosubtract(end+1,1) = 1;
                                           nocue_tosubtract(end,2) = cue_ar_no/250; %in seconds
                                           nocue_tosubtract(end,3) = a;
                                           nocue_tosubtract(end,4) = h;
                                end 
                            end
                        end
                    
                        
            end

            end
 toc
 
 
 %% Calculate total duration of the cue and no-cue period
 
 % cue
        cont_cue_duration = sum(period);
        cont_cue_duration = cont_cue_duration / 250; % in seconds
        cont_cue_duration = cont_cue_duration / 60; % in min

        indi_cues_duration = sum(sto_p) * cue_period;
        indi_cues_duration = indi_cues_duration / 250; % in seconds
        indi_cues_duration = indi_cues_duration / 60; % in min

        total_cueing_time = cont_cue_duration + indi_cues_duration;
        
        % Subtract arousals 
        sum_cue_tosubtract = sum(cue_tosubtract(:,2)); %in sec
        sum_cue_tosubtract = sum_cue_tosubtract/60; %in min
        
        total_cueing_time = total_cueing_time - sum_cue_tosubtract; %without arousals

  % no cue
        cont_nocue_duration = sum(period_nocue);
        cont_nocue_duration = cont_nocue_duration / 250; % in seconds
        cont_nocue_duration = cont_nocue_duration / 60; % in min
        
        total_break_time = cont_nocue_duration;
        
        % Subtract arousals
        sum_nocue_tosubtract = sum(nocue_tosubtract(:,2)); %in sec
        sum_nocue_tosubtract = sum_nocue_tosubtract/60; %in min
        
        total_break_time = total_break_time - sum_nocue_tosubtract; %without arousals

        
%% Save data
        ch_savePath = sprintf('%s/Stage%d',ch_savePath,stage);        
        ch_saveFile	= sprintf('%s_CuePeriodDuration%d.mat',mx_pFiles{pp},stage);       
        ch_saveFile = fullfile(ch_savePath,ch_saveFile);
        
        ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueingPeriod';

        fprintf('Saving file %s: ',ch_saveFile)
        tic
        save(ch_saveFile,'total_cueing_time','total_break_time','vt_cues');
        toc
        
        clear total_cueing_time total_break_time period_nocue vt_cues period sto_p epochduration vt_dHypno markers cues stageID timeID timeIDsamples correct_epochs correct_times cue_tosubtract nocue_tosubtract
 
end
