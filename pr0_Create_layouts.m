%% Create layouts for SP_SW detection

% choose your channels
vt_chEEG = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
              'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
              'P3','Pz','P4','P8', 'O1','O2'};
% Use a template from easycap 
cfg = [];
cfg.layout = 'easycapM1.mat'; %the easycapM1.mat file is available at 
                              %https://github.com/fieldtrip/fieldtrip.git 
                              %(got to template/layout/easycapM1.mat and 
                              %add the downloaded file to your path)
lay = ft_prepare_layout(cfg);

% Remove the channels that we dont have
remove = ismember(lay.label,vt_chEEG');
lay.pos = lay.pos(remove,:);
lay.height = lay.height(remove);
lay.label = lay.label(remove);
lay.width = lay.width(remove);

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
ft_layoutplot(cfglay);

save Layout_selected.mat lay %Layout will save as 'lay' in your current directory

clear('cfglay');
      