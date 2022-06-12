% Compute sleep spindles, slow oscillations and power spectrum density(PSD)

clc,clear

%% Define path
ch_hypnPath	= '/yourDirectory/Data/Sleep_scoring'; %sleep hypnograms folder 
ch_rootPath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/Preprocessed4spindles'; %preprocessed data folder
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/AllSpindles'; %output folder (needs to be created first)

%Define all participants
mx_pFiles	= [{'MRI_part7_sleep2'}]; %your participant of choice

%Define sleep stage
sleep_stage = [2,3]; % <----- remember to put the stage you want to analyse


        if sleep_stage == [2,3]
            stage = 23;         
        elseif sleep_stage == [2]
            stage = 2;
        elseif sleep_stage == [3]
            stage = 3;
        end

%% Define eeg variables

vt_chEEG    = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
               'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
               'P3','Pz','P4','P8','O1','O2'}; % all channels
vt_chEEGL   = {'FC3','C5','C3','C1','CP3','Fp1','F7','F3','T7','P7','P3','O1'}; % left hemisphere channels
vt_chEEGR   = {'FC4','C6','C4','C2','CP4','Fp2','F4','F8','T8','P4','P8','O2'}; % right hemisphere channels


%% Detection Settings
% Values based on Iber2007

nm_fSample  = 250;% sampling rate

vt_tfLims	= [0.1,35]; % general limits

vt_fPassEEG	= [0.3,35]; %not used
vt_fStopEEG	= [0.1,40]; %not used

vt_fPassSO	= [0.3,2]; %SO frequency
vt_fStopSO	= [0.1,4]; %bandpass filter
nm_troughThres  = -35; %threshold for the troughs

vt_fPassSp	= [11,16]; %spindles frequency
vt_fStopSp	= [9,18]; %bandpass filter

vt_timeSpindles	= [0.5,2];%duration of spindles
vt_timeFreqEvn	= [2,60]; %not used

nm_minNumOsc	= 5; % 5 oscialltions during spindle event
nm_scoreWind	= 30; % not used

nm_numFreq      = 256;
vt_freq         = linspace(vt_tfLims(1),vt_tfLims(2),nm_numFreq); % enerates linearly spaced vector (256 points between 0.1 and 35)

%% Create filters
fprintf('Design filters: ')
tic
ob_fSO	= fn_designIIRfilter(nm_fSample,vt_fPassSO,vt_fStopSO);
ob_fSP	= fn_designIIRfilter(nm_fSample,vt_fPassSp,vt_fStopSp);
toc

%% Define statistic variables (not used)
nm_maxStat  = 2;
nm_maxPerc	= 100*erf(nm_maxStat/sqrt(2));
nm_maxLossCh= 0.20; % maximum ratio of channels lost by trial
nm_numPerm  = 200;
nm_maxCumsum= 0.8;
nm_maxPercCh= 0.3;

%% Process data
for pp = 1:size(mx_pFiles,2) % loop for participant 

        st_events   = struct;
        ch_curFile	= sprintf('%s_SpindlesPreproc_CORRECT.mat',mx_pFiles{pp}); 
        ch_filename = fullfile(ch_rootPath,ch_curFile);
        
        nm_isfile   = exist(ch_filename,'file') == 2;
                               
        %% Load data
        fprintf('Loading file %s: ',ch_filename)
        tic
        st_dat	= load(ch_filename);      
        toc
        
        %% Load hypnogram
        
        fprintf('Loading data hypnogram: ')
        tic
        ch_hypfile  = sprintf('psgHypno-%s.mat',mx_pFiles{pp});
        st_hyp      = load(fullfile(ch_hypnPath,ch_hypfile));
        vt_dHypno   = single(st_hyp.dat(1,:));
        vt_dHypno   = interp1(st_hyp.timeEpoch,vt_dHypno,...
                    st_dat.time{1},'previous','extrap');
        toc
                
        %% Compute channel features
        st_events.so	= cell(size(vt_chEEG));
        st_events.sp	= cell(size(vt_chEEG));
        st_events.sPSD	= cell(size(vt_chEEG));
        st_events.mPSD	= cell(size(vt_chEEG));
        st_events.ePSD	= cell(size(vt_chEEG));
        st_events.fPSD	= vt_freq;
     
        for ch = 1:numel(vt_chEEG)
            %% Obtain time-frequency transform
            nm_curCh    = ismember(st_dat.label,vt_chEEG{ch});
            vt_signalCh	= st_dat.trial{1}(nm_curCh,:)';
            
            %% Filter in SO frequency band
            fprintf('Filtering in SO band for %s: ',vt_chEEG{ch})
            tic
            vt_signalSO	= fn_filterOffline(vt_signalCh,ob_fSO);
            toc
            
            %% Detect SO events
            fprintf('   ** Detect SO events: ')
            tic
            st_cnf              = struct;
            st_cnf.freqband     = [];
            st_cnf.fsampling	= st_dat.fsample; % samples
            st_cnf.threshold	= nm_troughThres; % trough threshold
            st_cnf.hypnogram    = vt_dHypno; % hypnogram
            st_cnf.stage        = sleep_stage; % sleep stage of interest
            st_cnf.minthresh	= [];
            st_cnf.toFilter     = [];
            
            st_events.so{ch}	= fn_detectsleepSO(vt_signalSO,st_cnf); %detects throughs
            toc
            
            clear st_cnf  
            
            %% Filter in the fast spindle band
            fprintf('Filtering in spindle band for %s: ',vt_chEEG{ch})
            tic
            vt_rmsFS	= single(fn_filterOffline(vt_signalCh,ob_fSP));% filters FS band using SP filter ob_fSP
            vt_rmsFS    = single(fn_rmstimeseries(vt_rmsFS,vt_timeSpindles(1)));% computes the root mean square timeseries of 
                                                                                 % signal vt_rmsFS using an sliding 
                                                                                 % window of vt_timeSpindles samples
            toc
            
            %% Detect spindle events            
            fprintf('	** Processing spindles: ')
            tic
            st_cnf              = struct;
            st_cnf.fsampling	= st_dat.fsample; %samples
            st_cnf.minnumosc	= nm_minNumOsc; % number of oscillations per spindle
            st_cnf.timebounds	= vt_timeSpindles; % spindle duration
            st_cnf.rawEEG       = vt_signalCh;
            st_cnf.freqband     = vt_fPassSp;
            st_cnf.method       = 'fixed';
            st_cnf.hypnogram    = vt_dHypno;
            st_cnf.stage        = sleep_stage;
            
            st_events.sp{ch}	= fn_detectsleepSpindles(vt_rmsFS,st_cnf);
            toc
            
            clear st_cnf vt_rmsFS
            
            %% Compute Ch PSD (power spectrum density)        
            fprintf('	** Processing psd: ')
            tic
            vt_hypno    = single(st_hyp.dat(1,:));
            vt_idSleep	= find(vt_hypno >= 1 & vt_hypno <= 5);% stage 1, 2, 3, REM (excludes wake)
            
            mx_psdStage = nan(numel(vt_freq),numel(vt_idSleep));% create samples by epoches matrix
            vt_psdStage = nan(1,numel(vt_idSleep));
            vt_psdEpoch = nan(1,numel(vt_idSleep));
            
            for ee = 1:numel(vt_idSleep)
                if vt_idSleep(ee) == numel(vt_hypno)
                    continue
                end
                
                nm_beg  = st_hyp.timeEpoch(vt_idSleep(ee)); %beginning
                nm_end  = st_hyp.timeEpoch(vt_idSleep(ee)+1); %ending
                
                vt_tPSD = st_dat.time{1} >= nm_beg & ...
                        st_dat.time{1} < nm_end; 
                    
                vt_sPSD = vt_signalCh(vt_tPSD);
                
                [vt_p,vt_f] = pwelch(vt_sPSD,[],[],[],nm_fSample);% Welch's power spectral density estimate
                vt_p        = 10*log10(vt_p);
                vt_psd      = interp1(vt_f,vt_p,vt_freq,'linear');
                
                mx_psdStage(:,ee)	= vt_psd(:);
                vt_psdStage(ee)     = vt_hypno(vt_idSleep(ee));
                vt_psdEpoch(ee)     = nm_beg;
            end
            
            % Check noisy stages
            vt_inPSD    = true(size(vt_psdStage));
            
            for ee = [1,2,3,5]
                vt_idStage	= find(vt_psdStage == ee);
                mx_psd      = mx_psdStage(:,vt_idStage);
                mx_meanPSD  = repmat(mean(mx_psd,2),1,numel(vt_idStage));
                mx_stdPSD   = repmat(std(mx_psd,[],2),1,numel(vt_idStage));
                vt_outPSD   = sum(mx_psd > mx_meanPSD + 3*(mx_stdPSD) | ...
                            mx_psd < mx_meanPSD - 3*(mx_stdPSD));
                vt_outPSD   = (vt_outPSD / numel(vt_freq)) > 1/4;
                
                vt_idStage  = vt_idStage(vt_outPSD);
                
                vt_inPSD(vt_idStage)	= false;
                clear vt_tPSD vt_tPSD vt_psd vt_p vt_f vt_idSleep
            end
            
            st_events.mPSD{ch}	= single(mx_psdStage(:,vt_inPSD));
            st_events.sPSD{ch}	= int8(vt_psdStage(:,vt_inPSD));
            st_events.ePSD{ch}	= single(vt_psdEpoch(:,vt_inPSD));
            clear mx_psdStage vt_psdStage vt_hypno vt_inPSD vt_outPSD
            clear vt_idStage mx_stdPSD mx_meanPSD mx_psd vt_signalCh
            toc
        
        
            %% Compute spindles density (spindles/min)
            nm_epochs = ismember(st_hyp.dat,sleep_stage); 
            nm_epochs = sum(nm_epochs); % nb of epochs in the stages of intered
            nm_min    = nm_epochs/2; % min spent in the stages of interest
            nm_events = size(st_events.sp{1,ch}); % number of spindles            
            nm_events = nm_events(1,1);
            spindles_density = nm_events/nm_min;% spindle density
            
           st_events.spDensity{ch} = spindles_density;% save spindle density for that channel
        
        end   
        
         %% Average spDensity for all channels, left and right hemisphere
        st_events.spAvDensityAcrossCH = mean(cell2mat(st_events.spDensity));
        vt_idEEGL	= ismember(st_dat.label,vt_chEEGL);
        vt_idEEGR	= ismember(st_dat.label,vt_chEEGR);
        st_events.spAvDensityAcrossLeftCH = mean(cell2mat(st_events.spDensity(vt_idEEGL))); 
        st_events.spAvDensityAcrossRightCH = mean(cell2mat(st_events.spDensity(vt_idEEGR))); 
        
        %% Saving file
               
        ch_savePath = sprintf('%s/Stage%d',ch_savePath,stage);        
        ch_saveFile	= sprintf('%s_Spindles_stage%d_allSpindlesAllCH.mat',mx_pFiles{pp},stage);       
        ch_saveFile = fullfile(ch_savePath,ch_saveFile);
        
        ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/AllSpindles'; %output folder

        fprintf('Saving file %s: ',ch_saveFile)
        tic
        save(ch_saveFile,'st_events','vt_chEEG');
        toc
        %clear st_events st_dStim st_comp st_dNoisy st_dTrial
        
   
end