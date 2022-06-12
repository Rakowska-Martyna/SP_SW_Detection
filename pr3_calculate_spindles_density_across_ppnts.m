% Calculate spindle density across participants

clc,clear

%% Define path
        
ch_savePath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/AllSpindles/AcrossPpnts';%output
ch_rootPath	= '/yourDirectory/Data/Sleep_EEG_analysis/Spindles/AllSpindles'; %data folder (pr2 output)

%Define all participants
 mx_pFiles	= [{'MRI_part7_sleep2'}];

nm_ppnt = numel(mx_pFiles);
SpindleResults = zeros(nm_ppnt,3); % create a ppnt by 3 matrix

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

for pp = 1:size(mx_pFiles,2) % loop for participant 

        st_events   = struct;
        ch_curFile	= sprintf('%s_Spindles_stage%d_allSpindlesAllCH.mat',mx_pFiles{pp},stage);
        ch_filename = fullfile(ch_rootPath,ch_curFile);
                               
        %% Load data
        fprintf('Loading ppnt %s from %s: \n',mx_pFiles{pp}, ch_filename)
        st_dat	= load(ch_filename);            
        
        %% Get spindle parameters
        
        avDensity = st_dat.st_events.spAvDensityAcrossCH;
        avLeftDensity = st_dat.st_events.spAvDensityAcrossLeftCH;
        avRightDensity = st_dat.st_events.spAvDensityAcrossRightCH;
        
        SpindlesResults(pp,1) = avDensity; % 1st column = average spindle density
        SpindlesResults(pp,2) = avLeftDensity; % 2nd column = left hemisphere spindle density
        SpindlesResults(pp,3) = avRightDensity; % 3rd column = right hemisphere spindle density
        
end

     ch_saveFile	= sprintf('SpindlesResults_stage%d.mat',stage);
     ch_saveFile = fullfile(ch_savePath,ch_saveFile);
         
     
    fprintf('Saving file %s: \n',ch_saveFile)
     save(ch_saveFile,'SpindlesResults');
