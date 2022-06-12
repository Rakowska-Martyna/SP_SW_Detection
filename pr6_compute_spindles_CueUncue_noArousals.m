%% Compute spindles density without arousals & during cue/no cue period (as in James 2014)
clc,clear

%% Define sleep stage of interest

sleep_stage = [2,3]; % <------ put the sleep stage you want to analyse

        if sleep_stage == [2,3]
            stage = 23;         
        elseif sleep_stage == [2]
            stage = 2;
        elseif sleep_stage == [3]
            stage = 3;
        end
        
%% Define path       
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueNocueSpindles'; %output directory
ch_rootPath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/noArousals'; %directory with idnetified spindles 
ch_rootPath = sprintf('%s/Stage%d',ch_rootPath,stage);  
ch_hypnPath	= '/yourDirectory/Data/Data/Sleep_scoring';%hypnogram direcory
ch_cuesPath = '/yourDirectory/Data/Sleep_EEG';%EEG directory
ch_cuePeriodPath = '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueingPeriod';%cue/no-cue period directory

%% Define participants
% Participnts from normal system (downsampled = 1) and from liveamp (not donwnsampled = 0)

downsampled = 1; % 1 = downsampled (500Hz -> 250Hz), 0 = not downsampled (250Hz)

if downsampled == 1
    mx_pFiles	= [{'MRI_part7_sleep2'}];
       
elseif downsampled == 0
    mx_pFiles	= [];
                
end

nm_ppnt = numel(mx_pFiles);
SpindleResults = zeros(nm_ppnt,3);

%% Define channels of interest
vt_chEEG    = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
               'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
               'P3','Pz','P4','P8','O1','O2'}; % all channels
vt_chEEGL   = {'FC3','C5','C3','C1','CP3','Fp1','F7','F3','T7','P7','P3','O1'}; % left hemisphere channels
vt_chEEGR   = {'FC4','C6','C4','C2','CP4','Fp2','F4','F8','T8','P4','P8','O2'}; % right hemisphere channels

cue_period = 875; % = 3.5 seconds * 250 Hz (sampling rate after downsampling) 
break_period = 5000; % 20 * 250Hz

%% Loop for participants
for pp = 1:size(mx_pFiles,2)                          
      
        %% Load data
        
        st_events   = struct;
        ch_curFile	= sprintf('%s_Spindles_stage%d_noArousals.mat',mx_pFiles{pp},stage);
        ch_filename = fullfile(ch_rootPath,ch_curFile);
        
        nm_isfile   = exist(ch_filename,'file') == 2;
        
        fprintf('Loading ppnt %s from %s: \n',mx_pFiles{pp}, ch_filename)
        st_dat	= load(ch_filename);            
        
        %% Load hypnogram       
        fprintf('Loading data hypnogram: ')
        
        tic
        ch_hypfile  = sprintf('psgHypno-%s.mat',mx_pFiles{pp});
        st_hyp      = load(fullfile(ch_hypnPath,ch_hypfile));
        vt_dHypno   = single(st_hyp.dat(1,:));
        vt_hArous   = st_hyp.arousals{1,1};
        toc
        
        %% Load cues (when the sounds were played)       
        fprintf('Loading cue markers: \n')
        
        ch_cuefile     = sprintf('%s.vmrk',mx_pFiles{pp}); % marker file      
        ch_cuefilename = fullfile(ch_cuesPath,ch_cuefile); 
        markers = importdata(ch_cuefilename,'',10000); 
        cues = markers(11:end,:); % the first few 'cues' isnt important  
        for k = 1:numel(cues)
            x = cues(k,1);
            y = strsplit(x{1,1},','); % vmrk file is in text file, need to split the string of characters... 
            vt_cues(k,1) = y(:,3); % ...and select only the 3rd 'word' (the actual sample)
        end
        
        vt_cues = str2double(vt_cues); 
        
        if downsampled == 1
        vt_cues = vt_cues/2; % need to divide by 2 for the downsampled data because the marker file was created in the original sampling rate
        end    
        
        %% Load cue period (duration of the cue / no cue period of the sleep stage of interest)
        
        ch_cuePeriodPathstage = sprintf('%s/Stage%d',ch_cuePeriodPath,stage);
        ch_curFile	= sprintf('%s_CuePeriodDuration%d.mat',mx_pFiles{pp},stage);
        ch_filename = fullfile(ch_cuePeriodPathstage,ch_curFile);
        cue_period_duaration = load(ch_filename);
        total_cueing_time = cue_period_duaration.total_cueing_time; % cue period
        total_break_time = cue_period_duaration.total_break_time; % no cue period
 
        %% Get spindles during cue/no cue periods
        
       fprintf('Getting spindles during cue/no cue periods \n')
        st_dat.st_dat.st_events.sp_cues = st_dat.st_dat.st_events.sp_noArousals; % create new variable for cue spindles that currently has all cue spindles excluding arousals
        st_dat.st_dat.st_events.sp_nocues = st_dat.st_dat.st_events.sp_noArousals;% create new variable for no-cue spindles that currently has all no-cue spindles excluding arousals
    
    for ch = 1:numel(vt_chEEG)
        fprintf('channel %d \n', ch)
        cue_sp{1,ch} = ones(size(st_dat.st_dat.st_events.sp_noArousals{1,ch})); % create 2 columns of 1s, size of all spindles
        nocue_sp{1,ch} = ones(size(st_dat.st_dat.st_events.sp_noArousals{1,ch})); % create 2 columns of 1s, size of all spindles
        events_sp = st_dat.st_dat.st_events.sp_noArousals{1,ch};
        nm_events = size(events_sp); % how many spindles without arousals
        nm_events = nm_events(1,1);
        nm_cues    = numel(vt_cues); % how many cues
        for i = 1:nm_events % for each spindle
        % cues
            for a = 1:nm_cues % for each cue
                x = events_sp(i,1):events_sp(i,2); % all samples within a spindle
                for j = 1:numel(x) % for each sample within a spindle
                if any(x(j) >= vt_cues(a,1) && x(j) <= (vt_cues(a,1) + cue_period)) % if any sample of a chosen spindle is between the start of a cue and [cue + cue_period]                  
                   cue_sp{1,ch}(i,:) = 0; % this spindle is a cue spindle, cue_sp is 0 for that spindle
                % no cues
                elseif a ~= nm_cues
                    if any(x(j) >= (vt_cues(a,1) + cue_period) && x(j) <= (vt_cues(a,1) + break_period) && (vt_cues(a+1,1)-vt_cues(a,1)) >= break_period && (vt_cues(a+1,1)-vt_cues(a,1)) <= (cue_period + break_period)) % if any sample of a chosen spindle is between [cue + cue_period] and [cue + break_period]
                        nocue_sp{1,ch}(i,:) = 0;  % this spindle is a no-cue spindle, nocue_sp is 0 for that spindle
                    end
                    
                end
                end
            end
        end
    end
  
    
           %% Get SOs during cue/no cue periods
        
      fprintf('Getting SOs during cue/no cue periods \n')
        st_dat.st_dat.st_events.so_cues = st_dat.st_dat.st_events.so_noArousals; 
        st_dat.st_dat.st_events.so_nocues = st_dat.st_dat.st_events.so_noArousals;
    
        for ch = 1:numel(vt_chEEG)
        fprintf('channel %d \n', ch)
        cue_so{1,ch} = ones(size(st_dat.st_dat.st_events.so_noArousals{1,ch}));
        nocue_so{1,ch} = ones(size(st_dat.st_dat.st_events.so_noArousals{1,ch}));
        events_so = st_dat.st_dat.st_events.so_noArousals{1,ch};
        nm_so_events = size(events_so);
        nm_so_events = nm_so_events(1,1);
        nm_cues    = numel(vt_cues);
        for i = 1:nm_so_events
        % cues
            for a = 1:nm_cues
                if any(events_so(i) >= vt_cues(a,1) && events_so(i) <= (vt_cues(a,1) + cue_period))                 
                   cue_so{1,ch}(i,:) = 0; % if there is a SO during cue, cue_sp is 0
        % no cues
                elseif a ~= nm_cues
                    if any(events_so(i) >= (vt_cues(a,1) + cue_period) && events_so(i) <= (vt_cues(a,1) + break_period) && (vt_cues(a+1,1)-vt_cues(a,1)) >= break_period && (vt_cues(a+1,1)-vt_cues(a,1)) <= (cue_period + break_period))
                      nocue_so{1,ch}(i,:) = 0; 
                    end
                end
            end
        end
        end
      
    %% end
    
    for ch = 1:numel(vt_chEEG) 
        
        % SPINDLES
        
        % if a spindle is both in cue and no-cue period (falls in between) then exclude it
        for f=1:length(cue_sp{1,ch}(:,1))
            if sum(cue_sp{1,ch}(f,:) == nocue_sp{1,ch}(f,:))>1
                cue_sp{1,ch}(f,:) = 1;
                nocue_sp{1,ch}(f,:) = 1;
            end
        end
        

        st_dat.st_dat.st_events.sp_cues{1,ch} = st_dat.st_dat.st_events.sp_cues{1,ch}(~cue_sp{1,ch}(:,1),:);       
        st_dat.st_dat.st_events.sp_nocues{1,ch} = st_dat.st_dat.st_events.sp_nocues{1,ch}(~nocue_sp{1,ch}(:,1),:);       
        
        % SOs
        % if SO is both in cue and no-cue period (falls in between) then exclude it
          for f=1:length(cue_so{1,ch}(:,1))
            if sum(cue_so{1,ch}(f,:) == nocue_so{1,ch}(f,:))>1
                cue_so{1,ch}(f,:) = 1;
                nocue_so{1,ch}(f,:) = 1;
            end
           end
        
        st_dat.st_dat.st_events.so_cues{1,ch} = st_dat.st_dat.st_events.so_cues{1,ch}(~cue_so{1,ch}(:,1),:);       
        st_dat.st_dat.st_events.so_nocues{1,ch} = st_dat.st_dat.st_events.so_nocues{1,ch}(~nocue_so{1,ch}(:,1),:);       
        
        %% Compute spindles density   
        % spindles/duration of cue/nocue period during the stage of interest 
        % =total_cueing_time & total break time
         
        % compute number of events
        % cued
        nm_events = size(st_dat.st_dat.st_events.sp_cues{1,ch});           
        nm_events = nm_events(1,1);
        % uncued
        nm_events_nocues = size(st_dat.st_dat.st_events.sp_nocues{1,ch});           
        nm_events_nocues = nm_events_nocues(1,1);
              
        % compute spindles density
        % cued
        spindles_density = nm_events/total_cueing_time;            
        st_dat.st_dat.st_events.spDensity_cue{ch} = spindles_density;
        % uncued
        spindles_density_nocues = nm_events_nocues/total_break_time;            
        st_dat.st_dat.st_events.spDensity_nocue{ch} = spindles_density_nocues;
        
        %% Compute SOs     
        
        % compute number of events
        % cued
        nm_so_events = size(st_dat.st_dat.st_events.so_cues{1,ch});           
        nm_so_events = nm_so_events(1,1);
        % uncued
        nm_so_events_nocues = size(st_dat.st_dat.st_events.so_nocues{1,ch});           
        nm_so_events_nocues = nm_so_events_nocues(1,1);
                
        % compute SOs density
        % cued
        so_density = nm_so_events/total_cueing_time;            
        st_dat.st_dat.st_events.soDensity_cue{ch} = so_density;
        % uncued
        so_density_nocues = nm_so_events_nocues/total_break_time;            
        st_dat.st_dat.st_events.soDensity_nocue{ch} = so_density_nocues;
         
    end
    
        %% Summary stats = spidnles
        %Summary stats for cue period  
        st_dat.st_dat.st_events.spAvDensityAcrossCH_cue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_cue));
        vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL); % left electrodes
        vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR); % right electrodes
        st_dat.st_dat.st_events.spAvDensityAcrossLeftCH_cue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_cue(vt_idEEGL))); 
        st_dat.st_dat.st_events.spAvDensityAcrossRightCH_cue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_cue(vt_idEEGR))); 
        
        % Summary stats for NOcue period
         st_dat.st_dat.st_events.spAvDensityAcrossCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_nocue));
         vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL); % left electrodes
         vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR); % right electrodes
         st_dat.st_dat.st_events.spAvDensityAcrossLeftCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_nocue(vt_idEEGL))); 
         st_dat.st_dat.st_events.spAvDensityAcrossRightCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.spDensity_nocue(vt_idEEGR))); 
        
        %% Summary stats = SOs
        %Summary stats for cue period        
        st_dat.st_dat.st_events.soAvDensityAcrossCH_cue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_cue));
        vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL); % left electrodes
        vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR); % right electrodes
        st_dat.st_dat.st_events.soAvDensityAcrossLeftCH_cue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_cue(vt_idEEGL))); 
        st_dat.st_dat.st_events.soAvDensityAcrossRightCH_cue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_cue(vt_idEEGR))); 
          
        % Summary stats for NOcue period
         st_dat.st_dat.st_events.soAvDensityAcrossCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_nocue));
         vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL); % left electrodes
         vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR); % right electrodes
         st_dat.st_dat.st_events.soAvDensityAcrossLeftCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_nocue(vt_idEEGL))); 
         st_dat.st_dat.st_events.soAvDensityAcrossRightCH_nocue = mean(cell2mat(st_dat.st_dat.st_events.soDensity_nocue(vt_idEEGR))); 
                
        %% Save data
        
        ch_savePath = sprintf('%s/Stage%d',ch_savePath,stage);        
        ch_saveFile	= sprintf('%s_Spindles_stage%d_CueNoCue_NEW.mat',mx_pFiles{pp},stage);       
        ch_saveFile = fullfile(ch_savePath,ch_saveFile);
        
        ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueNocueSpindles';

        fprintf('Saving file %s: ',ch_saveFile)
        tic
        save(ch_saveFile,'st_dat','vt_chEEG');
        toc
        
        clear vt_cues
        
end

    
    