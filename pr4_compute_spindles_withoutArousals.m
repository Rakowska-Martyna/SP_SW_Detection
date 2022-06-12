% Spindles density without arousals
clc,clear

% (1) Multiplies arousals (start & end) by sampling rate (250Hz for all participants)
% (2) Removes spindles that fall within an arousal 
%     (i.e. if any part of a spindle falls within any arousal that spindle is removed)
% (3) Calculates spindle density by dividing the number of remaining spindles
%     by the total time in the chosen sleep stage, minus the total duration of
%     all arousals in the chosen sleep stage
% (4) Calculates SOs density in the same manner

%% Define path
        
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/noArousals'; %output folder
ch_rootPath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/AllSpindles'; %data folder 
ch_hypnPath	= '/yourDirectory/Data/Sleep_scoring';

% Define ppnts
mx_pFiles	= [{'MRI_part7_sleep2'}];

nm_ppnt = numel(mx_pFiles);
SpindleResults = zeros(nm_ppnt,3);

% Define sleep stage
sleep_stage = [2,3]; % <------ put the stage you want to analyse
        if sleep_stage == [2,3]
            stage = 23;         
        elseif sleep_stage == [2]
            stage = 2;
        elseif sleep_stage == [3]
            stage = 3;
        end
        
ch_rootPath = sprintf('%s/Stage%d',ch_rootPath,stage);
          
%% Define eeg variables
vt_chEEG    = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
               'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
               'P3','Pz','P4','P8','O1','O2'}; % all channels
vt_chEEGL   = {'FC3','C5','C3','C1','CP3','Fp1','F7','F3','T7','P7','P3','O1'}; % left hemisphere channels
vt_chEEGR   = {'FC4','C6','C4','C2','CP4','Fp2','F4','F8','T8','P4','P8','O2'}; % right hemisphere channels


remove_sp = cell(1,(numel(vt_chEEG)));
remove_so = cell(1,(numel(vt_chEEG)));

%% Remove spindles that fall within arousals

% Loop for participant 
for pp = 1:size(mx_pFiles,2) 
    
        %% Load data        
        st_events   = struct;
        ch_curFile	= sprintf('%s_Spindles_stage%d_allSpindlesAllCH.mat',mx_pFiles{pp},stage);
        ch_filename = fullfile(ch_rootPath,ch_curFile);
        
        fprintf('Loading ppnt %s from %s: \n',mx_pFiles{pp}, ch_filename)
        st_dat	= load(ch_filename);            
        
         %% Load hypnogram        
        fprintf('Loading data hypnogram \n')
        ch_hypfile  = sprintf('psgHypno-%s.mat',mx_pFiles{pp});
        st_hyp      = load(fullfile(ch_hypnPath,ch_hypfile));
        vt_dHypno   = single(st_hyp.dat(1,:));
        vt_hArous   = st_hyp.arousals{1,1};
         
        vt_hArous = vt_hArous * 250; % convert seconds to samples
        
        %% Remove spindles & SO during arousal        
        arousal.start = vt_hArous(:,1); 
        arousal.end = vt_hArous(:,2);
    
    % SPINDLES
    fprintf('Removing spindles during arousals \n')
    st_dat.st_events.sp_noArousals = st_dat.st_events.sp;
    for ch = 1:numel(vt_chEEG) % for each channel
        remove_sp{1,ch} = zeros(size(st_dat.st_events.sp{1,ch})); % create an empty array to mark spindles for removal
        events_sp = st_dat.st_events.sp{1,ch}; % all the spindles detected before
        nm_events = size(events_sp);
        nm_events = nm_events(1,1);
        nm_arousal = numel(arousal.start);
        for i = 1:nm_events % for all the spindles detected before
            for k = 1:nm_arousal % and all the arousals               
                x = events_sp(i,1):events_sp(i,2); % for all the samples within the spindle
                for j = 1:numel(x)
                if any(x(j) >= arousal.start(k,1) && x(j) <= arousal.end(k,1)) % find if any sample within a spindles is within an arousal
                %fprintf 'spindle during arousal \n';
                remove_sp{1,ch}(i,:) = 1; % if there is, mark that spindle for removal
                end             
                end              
            end
        end
    end
    
    % SLOW OSCILLATIONS
    fprintf('Removing SOs during arousals \n')
    st_dat.st_events.so_noArousals = st_dat.st_events.so;
    for ch = 1:numel(vt_chEEG) % for each channel
        remove_so{1,ch} = zeros(size(st_dat.st_events.so{1,ch})); % create an empty array to mark slow oscilations for removal
        events_so = st_dat.st_events.so{1,ch}; % all the slow oscilations detected before
        nm_so_events = size(events_so);
        nm_so_events = nm_so_events(1,1);
        nm_arousal = numel(arousal.start);
        for i = 1:nm_so_events % for all the slow oscillations detected before
            for k = 1:nm_arousal % and all the arousals               
                if events_so(i) >= arousal.start(k,1) && events_so(i) <= arousal.end(k,1) % find if the SO is within an arousal
                %fprintf 'SO during arousal \n';
                remove_so{1,ch}(i,:) = 1; % if there is, mark that SO for removal
                end                          
            end
        end
     end
            
    for ch = 1:numel(vt_chEEG)

        %SPINDLES
        st_dat.st_events.sp_noArousals{1,ch} = st_dat.st_events.sp_noArousals{1,ch}(~remove_sp{1,ch}(:,1),:);
         
        % SLOW OSCILLATIONS
        st_dat.st_events.so_noArousals{1,ch} = st_dat.st_events.so_noArousals{1,ch}(~remove_so{1,ch}(:,1),:);
        
        %%% Comput arousals time that needs to be subtracted from total time %%% 
        remove_epochs = zeros(length(st_hyp.arousals{1,1}),1);
        remove_sec = zeros(length(st_hyp.arousals{1,1}),1);
        for hyp = 1:length(st_hyp.arousals{1,1})
            length_ar = st_hyp.arousals{1,1}(hyp,2) - st_hyp.arousals{1,1}(hyp,1);
            epoch = (st_hyp.arousals{1,1}(hyp,1)/30)+1;
            epoch = floor(epoch); % round down to get the epoch number where it starts
            if sum(st_hyp.dat(epoch) == sleep_stage) == 1
                remove_epochs(hyp,1) = 1;
                remove_sec(hyp,1) = length_ar;
            end
            
        end
        
        sum_remove_epochs = sum(remove_epochs); 
        sum_remove_length = sum(remove_sec);
        sum_remove_epochs_min = sum_remove_epochs/2;
        sum_remove_length_min = sum_remove_length/60;
  
        %%Compute spindles density (spindles/min)
            nm_epochs = ismember(st_hyp.dat,sleep_stage); 
            nm_epochs = sum(nm_epochs); % nb of epochs in the stages of intered
            nm_min    = nm_epochs/2;
            nm_events = size(st_dat.st_events.sp_noArousals{1,ch});           
            nm_events = nm_events(1,1);
            spindles_density = nm_events/nm_min;
            spindles_density = nm_events/(nm_min-sum_remove_length_min);
            
            st_dat.st_events.spDensity_noArousals{ch} = spindles_density;
            
       %%Compute sos density (spindles/min)
            nm_epochs = ismember(st_hyp.dat,sleep_stage); 
            nm_epochs = sum(nm_epochs); % nb of epochs in the stages of intered
            nm_min    = nm_epochs/2;
            nm_so_events = size(st_dat.st_events.so_noArousals{1,ch});           
            nm_so_events = nm_so_events(1,1);
            so_density = nm_so_events/nm_min;
            so_density = nm_so_events/(nm_min-sum_remove_length_min);
            
            st_dat.st_events.soDensity_noArousals{ch} = so_density;
    end
        
        %% Summary stats
        % SPINDELS
        st_dat.st_events.spAvDensityAcrossCH_noArousals = mean(cell2mat(st_dat.st_events.spDensity_noArousals));
        vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL);
        vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR);
        st_dat.st_events.spAvDensityAcrossLeftCH_noArousals = mean(cell2mat(st_dat.st_events.spDensity_noArousals(vt_idEEGL))); 
        st_dat.st_events.spAvDensityAcrossRightCH_noArousals = mean(cell2mat(st_dat.st_events.spDensity_noArousals(vt_idEEGR))); 
        
        % SOs
        st_dat.st_events.soAvDensityAcrossCH_noArousals = mean(cell2mat(st_dat.st_events.soDensity_noArousals));
        vt_idEEGL	= ismember(st_dat.vt_chEEG,vt_chEEGL);
        vt_idEEGR	= ismember(st_dat.vt_chEEG,vt_chEEGR);
        st_dat.st_events.soAvDensityAcrossLeftCH_noArousals = mean(cell2mat(st_dat.st_events.soDensity_noArousals(vt_idEEGL))); 
        st_dat.st_events.soAvDensityAcrossRightCH_noArousals = mean(cell2mat(st_dat.st_events.soDensity_noArousals(vt_idEEGR))); 
        
        %% Save data
        ch_savePath = sprintf('%s/Stage%d',ch_savePath,stage);        
        ch_saveFile	= sprintf('%s_Spindles_stage%d_noArousals.mat',mx_pFiles{pp},stage);       
        ch_saveFile = fullfile(ch_savePath,ch_saveFile);
        
        ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/noArousals';

        fprintf('Saving file %s: ',ch_saveFile)
        tic
        save(ch_saveFile,'st_dat','vt_chEEG');
        toc
        
end  