% Preprocess EEG data for spindle analysis

clc,clear

%% Define path
ch_hypnPath	= 'yourDirectory/Data/Sleep_scoring'; %sleep hypnograms folder 
ch_rootPath	= 'yourDirectory/Data/Sleep_EEG'; %raw EEG data folder 
ch_savePath	= 'yourDirectory/Data/Sleep_EEG_analysis/Spindles/Preprocessed4spindles'; %output folder (needs to be created first)
% Load layout
load('/yourDirectory/Layouts/Layout_all.mat');
% Add Miguel's functions
addpath('/yourDirectory/https://github.com/mnavarretem/myToolbox.git');

%%  Choose data/system to analyse 
% Since half of the data was collected using liveAmp (sampling rate = 250Hz) and half using
% BrainAmp (the 'normal system', sampling rate = 500Hz) we need to
% downsample the latter and hence analyse the data separately

liveamp = 1; % 0 = normal system (500Hz, need to downsample), 1 = LiveAmp (250Hz)

if liveamp == 0
    mx_pFiles	= [{'MRI_part3_sleep'},{'MRI_part15_sleep1'}]; 
           
elseif liveamp == 1
    mx_pFiles	= [{'MRI_part7_sleep2'}];% a participant of your choice
 
end

%% Select particiapnts with bad ref channel(s) (i.e. channels impossible to score) 
% - bad channel interpolated before, now need to replace it here

int_ref_Path = 'yourDirectory/Sleep_EEG/Interpolated_channels' ; %interpolated channels folder
bad_ref_ppnt = {'MRI_part8_sleep.eeg','MRI_part16_sleep1.eeg'}; 
bad_ref_ppnt_filename = fullfile(ch_rootPath,bad_ref_ppnt);
bad_channel = 'TP10'; %interpolated channel

%% Define your EEG variables

vt_chEEG    = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
               'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
               'P3','Pz','P4','P8','O1','O2'}; % all channels
vt_chEEGL   = {'FC3','C5','C3','C1','CP3','Fp1','F7','F3','T7','P7','P3','O1'}; % left hemisphere channels
vt_chEEGR   = {'FC4','C6','C4','C2','CP4','Fp2','F4','F8','T8','P4','P8','O2'}; % right hemisphere channels

%% Define signal processing variables
nm_fSample  = 250; % because you will filter after downsampling so everything will already be 250Hz
vt_fPassEEG  = [0.3,35]; % pass band (accept anything between 0.3 and 35)
vt_fStopEEG  = [0.1,45]; % stop band (reject anything below 0.1 and above 45)

%% Create filters

fprintf('Design filters: ')
ob_fEEG	= fn_designIIRfilter(nm_fSample,vt_fPassEEG,vt_fStopEEG);

%% Preprocessing loop
for pp = 1:size(mx_pFiles,2)
      
        %% Select file       
        ch_file     = sprintf('%s.eeg',mx_pFiles{pp});
        ch_filename = fullfile(ch_rootPath,ch_file);
        [~,ch_file]	= fileparts(mx_pFiles{pp});

        %% Load data
        fprintf('Processing ppnt %s: \n',mx_pFiles{pp})
        fprintf('Loading file %s: \n',ch_filename)
        st_cfg                  = struct;
        st_cfg.dataset          = ch_filename;
        st_cfg.channel          = horzcat(vt_chEEG,vt_chEOG,vt_chEMG,vt_chREF);
        st_cfg.feedback         = 'no';
        st_cfg.trackcallinfo	= 'no';
tic
        st_dat	= ft_preprocessing(st_cfg); % reads EEG data to the memory 

        [~,vt_idDat]	= ismember(st_dat.label,st_cfg.channel);
        [~,vt_idDat]    = sort(vt_idDat);

        % Order st_dat.trial and st_dat.label in vt_idDat
        st_dat.trial{1} = st_dat.trial{1}(vt_idDat,:);
        st_dat.label	= st_dat.label(vt_idDat);
toc
       
  %% Replace bad ref with interpolated ref for chosen participants       
  
  % if the current ppnt is the one with bad red, find where the bad ref is
   if ismember(ch_filename,bad_ref_ppnt_filename)
       intch=[];
        for x = 1:length(st_dat.label)
            if ismember(st_dat.label(x,1),bad_channel)
                intch = x; % st_dat.label(x,1) is the bad channel = interpolated electrode
                cd(int_ref_Path)
                fprintf('Using interpolated ref')
            end
        end
     
   % go to the folder with the interpolated ref and find the correct ppnt
       cd(int_ref_Path)
       
       listing = dir;
    for K   = 1:length(listing)
        fname = listing(K).name;
        name  = mx_pFiles{pp}(1:15);
        if length(fname) > 3
            if strcmp(fname(1:15), name); break; end
        end
    end 
    
    % load the interpolated ref and replace the bad ref with the interpolated one
       load(fname)
       st_dat.trial{1,1}(intch,:) = st_dat_badchannel;
   end
     
        %% Down-sample 500Hz data to 250Hz 
        fprintf('Downsampling \n')
        if liveamp == 0
            st_datDS = st_dat;
            st_datDS.fsample  = 250; 
            n = length(st_dat.trial{1});        
            st_datDS.trial = cellfun(@(x) x(:,1:2:n),st_dat.trial,'UniformOutput',false);
            st_datDS.time = cellfun(@(x) x(:,1:2:n),st_dat.time,'UniformOutput',false);
            st_dat = st_datDS;
        end
        
        %% Filter data (EEG & REF)
        % Define ch/filter type
        vt_id                       = ismember(st_dat.label,vt_chEEG);
        st_dat.hdr.chantype(vt_id)	= {'eeg'};

        % Filter EEG data with the filter you crearted       
        vt_id  = ismember(st_dat.label,vt_chEEG) | ismember(st_dat.label,vt_chREF);
        fprintf('Filtering EEG: \n') 
        tic
        st_dat.trial{1}(vt_id,:)	= fn_filterOffline(...
                                    st_dat.trial{1}(vt_id,:)',ob_fEEG)';
        toc

         st_dat.label = st_dat.label(:); % Transpose the label to a column vector, otherwise you get an error during the channel repair
    
         %% Prepare layout
        % Use a template from easycap 
        cfg = [];
        cfg.layout = 'easycapM1.mat';
        lay = ft_prepare_layout(cfg);

        % Visualise easycapM1
        cfg = [];
        cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
        %ft_layoutplot(cfg);

        % Remove the channels that we dont have
        remove = ismember(lay.label,vt_chEEG');

        lay.pos = lay.pos(remove,:);
        lay.height = lay.height(remove);
        lay.label = lay.label(remove);
        lay.width = lay.width(remove);

        % Visualise our final layout
        cfg = [];
        cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
        %ft_layoutplot(cfg);
 
        %save LayoutSpindles_Exp1.mat lay
        %%% Load the channel information and prepare the layout
         elec.label       = lay.label;   % Wherever you have your label info
         elec.elecpos     = lay.pos;     % A matrix with the electrode locations
         elec.chanpos     = lay.pos;
         dat.elec         = elec;
         cfglay           = [];
         cfglay.layout    = lay;
         cfglay.elec      = elec;
         cfglay.rotate    = 90;          % The origin is different in Fieldtrip, so we have to rotate it
         cfglay.skipscale = 'yes';
         cfglay.skipcomnt = 'yes';
         layout           = ft_prepare_layout(cfglay, dat);
         ft_layoutplot(cfglay); % visualises final layout

         clear('cfglay');
         
         %% Exclude & interpolate bad channels

        %%% Preparing neighbours for channel interpolation and repair
        cfg          = [];
        cfg.method   = 'triangulation';  % calculates a triangulation based on a two-dimenstional projection of the sensor position (the layout)
        cfg.layout   = layout;           % the layout prepared above
        cfg.channel  = 'eeg';
        cfg.feedback = 'no';             % normally 'yes' - you have to check that the layout is correct
        neighbours   = ft_prepare_neighbours(cfg);

        %%% Summary Mode
        cfg             = [];
        cfg.method      = 'summary';
        cfg.neighbours  = neighbours;                % The neighbours prepared above; this is necessary for channel repair to work.
        cfg.keepchannel = 'repair';                 % Repair any rejected channels using ft_channelrepair
        cfg.layout      = layout;  
        cfg.elec        = dat.elec;
        cfg.trials      = 'all'; %trial_list;
        st_dat          = ft_rejectvisual(cfg, st_dat); % Browse channels. Make sure to click 'quit' at the end, otherwise it doesn't save them.

        % Just visualise channels (same y axis for all)
        for ff = 1:length(vt_chEEG)
        cfg = [];
        cfg.ylim = [-5000 5000];
        cfg.channel = vt_chEEG(ff);
        figure; ft_singleplotER(cfg,st_dat);
        end
        
        %%% Channel Mode - visualise and remove bad
        % Gives an overview of every trial for one channel at a time. 
        % Allows you to reject channels/electrodes that might have fallen off
        cfg             = [];
        cfg.method      = 'channel';
        cfg.neighbours  = neighbours;                % The neighbours prepared above; this is necessary for channel repair to work.
        cfg.keepchannel = 'repair';                 % Repair any rejected channels using ft_channelrepair
        cfg.layout      = layout;  
        cfg.elec        = dat.elec;
        cfg.trials      = 'all'; %trial_list;
        st_dat          = ft_rejectvisual(cfg, st_dat); % Browse channels. Make sure to click 'quit' at the end, otherwise it doesn't save them.
        
        %% Re-reference
        fprintf('Rereferencing \n')
        vt_refId     = ismember(st_dat.label,vt_chREF);
        vt_refAvg	= st_dat.trial{1}(vt_refId,:);
        vt_refAvg   = mean(vt_refAvg);
        vt_isEEG    = ismember(st_dat.label,vt_chEEG);
        mx_data     = st_dat.trial{1}(vt_isEEG,:);
        mx_data     = mx_data - repmat(vt_refAvg,size(mx_data,1),1);

        st_dat.trial{1}(vt_isEEG,:)	= mx_data;
        
        st_dat.trial{1}	= st_dat.trial{1}(~vt_refId,:); 
        st_dat.label 	= st_dat.label(~vt_refId);
        st_cfg.channel  = st_cfg.channel(~vt_refId);
        st_dat.hdr.chantype	= st_dat.hdr.chantype(~vt_refId);
        
        st_dat.time{1}  = single(st_dat.time{1});
        st_dat.trial{1} = single(st_dat.trial{1});

        %% Save clean file

        ch_saveFile	= sprintf('%s_SpindlesPreproc_CORRECT.mat',mx_pFiles{pp});
        ch_saveFile = fullfile(ch_savePath,ch_saveFile);

        fprintf('Saving file %s: \n',ch_saveFile)
        
        save(ch_saveFile,'-struct','st_dat');
       
        clear st_dat %mx_trials vt_clicks
    
end
