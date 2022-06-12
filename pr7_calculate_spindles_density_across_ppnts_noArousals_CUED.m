% spindles density across ppnts
clc,clear

%% Define sleep stage

sleep_stage = [2,3]; % <------ put the stage you want to analyse

        if sleep_stage == [2,3]
            stage = 23;         
        elseif sleep_stage == [2]
            stage = 2;
        elseif sleep_stage == [3]
            stage = 3;
        end
        
%% Define paths
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueNocueSpindles/AcrossPpnts';
ch_rootPath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/CueNocueSpindles';
ch_rootPath = sprintf('%s/Stage%d',ch_rootPath,stage);  

%% Define participants 

mx_pFiles	= [{'MRI_part7_sleep2'}];
nm_ppnt = numel(mx_pFiles);

CueSpindleResults = zeros(nm_ppnt,3);
NoCueSpindleResults = zeros(nm_ppnt,3);
CueSOsResults = zeros(nm_ppnt,3);
NoCueSOsResults = zeros(nm_ppnt,3);

%% Loop for participant
for pp = 1:size(mx_pFiles,2)                              
        
        %% Load data
        
        st_events   = struct;
        ch_curFile	= sprintf('%s_Spindles_stage%d_CueNoCue_NEW.mat',mx_pFiles{pp},stage);        
        ch_filename = fullfile(ch_rootPath,ch_curFile);        
        nm_isfile   = exist(ch_filename,'file') == 2;
        
        fprintf('Loading ppnt %s from %s: \n',mx_pFiles{pp}, ch_filename)
        st_dat	= load(ch_filename);            
        
        %% Get spindle parameters
        
        % Cue period
        avDensity = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossCH_cue; % average density (all motor channels)
        avLeftDensity = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossLeftCH_cue; % left channels density
        avRightDensity = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossRightCH_cue; % right channels density)
        
        CueSpindlesResults(pp,1) = avDensity;
        CueSpindlesResults(pp,2) = avLeftDensity;
        CueSpindlesResults(pp,3) = avRightDensity;
        
        % No cue period
        avDensity_noCue = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossCH_nocue; % average density (all motor channels)
        avLeftDensity_noCue = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossLeftCH_nocue; % left channels density
        avRightDensity_noCue = st_dat.st_dat.st_dat.st_events.spAvDensityAcrossRightCH_nocue; % right channels density)
        
        NoCueSpindleResults(pp,1) = avDensity_noCue;
        NoCueSpindleResults(pp,2) = avLeftDensity_noCue;
        NoCueSpindleResults(pp,3) = avRightDensity_noCue;
 
        %% Get SOs parameters
        
        % Cue period
        av_so_Density = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossCH_cue; % average density (all motor channels)
        avLeft_so_Density = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossLeftCH_cue; % left channels density
        avRight_so_Density = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossRightCH_cue; % right channels density)
        
        CueSOsResults(pp,1) = av_so_Density;
        CueSOsResults(pp,2) = avLeft_so_Density;
        CueSOsResults(pp,3) = avRight_so_Density;
        
        % No cue period
        av_so_Density_noCue = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossCH_nocue; % average density (all motor channels)
        avLeft_so_Density_noCue = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossLeftCH_nocue; % left channels density
        avRight_so_Density_noCue = st_dat.st_dat.st_dat.st_events.soAvDensityAcrossRightCH_nocue; % right channels density)
        
        NoCueSOsResults(pp,1) = av_so_Density_noCue;
        NoCueSOsResults(pp,2) = avLeft_so_Density_noCue;
        NoCueSOsResults(pp,3) = avRight_so_Density_noCue;
end

%% Save sumary data across participants

    ch_saveFile	= sprintf('SpindlesResults_stage%d_CueNoCue_NEW.mat',stage);
    ch_saveFile = fullfile(ch_savePath,ch_saveFile);
        
    fprintf('Saving file %s: \n',ch_saveFile)
    save(ch_saveFile,'CueSpindlesResults','NoCueSpindleResults','CueSOsResults','NoCueSOsResults');
