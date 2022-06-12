 %% Combine motor and other electrodes
  
 %Create topoplots from copied and pasted data
 
  %cueN2     = [4.28	5.06	4.49	4.23	4.53	3.88	4.21	4.42	4.29	4.48	4.59	4.37	4.79	4.72	4.86	4.67	4.48	4.03	4.00	3.73	4.27	4.37	4.43	4.37	3.85	4.25	3.33] ;
  %cueN3     = [3.70	4.20	4.21	3.32	3.69	3.28	3.76	3.54	3.92	3.87	3.98	3.62	4.13	4.19	4.14	4.16	3.71	3.28	3.49	2.97	3.30	3.93	4.00	3.64	3.12	3.55	2.97];
  cueN23    = [4.54	4.39	4.09	3.98	3.70	3.59	4.27	4.03	4.85	4.56	5.20	4.91	5.22	5.34	4.96	4.94	4.64	4.24	4.01	4.06	3.59	3.98	4.17	3.91	3.69	3.57	3.36];
  
  %nocueN2   = [3.54	4.22	3.67	3.42	3.79	3.12	3.49	3.50	3.67	3.81	3.83	3.48	3.87	4.02	3.96	3.92	3.68	3.38	3.27	3.07	3.52	3.72	3.77	3.46	3.31	3.36	2.75];
  %nocueN3   = [3.33	4.04	3.90	2.96	3.47	3.04	3.45	3.29	3.61	3.58	3.55	3.14	3.70	3.84	3.66	3.80	3.32	3.07	3.29	2.64	2.93	3.60	3.61	3.21	2.85	3.15	2.55];
  nocueN23  = [4.05	3.97	3.69	3.73	3.33	3.16	3.94	3.68	4.38	4.06	4.59	4.29	4.53	4.76	4.38	4.33	3.99	3.58	3.67	3.33	3.11	3.63	3.93	3.53	3.24	3.24	3.09];
  
  %allN2     = [3.91	4.64	4.08	3.83	4.16	3.50	3.85	3.96	3.98	4.15	4.21	3.93	4.33	4.37	4.41	4.30	4.08	3.70	3.64	3.40	3.90	4.04	4.10	3.91	3.58	3.80	3.04];
  %allN3     = [3.52	4.12	4.06	3.14	3.58	3.16	3.61	3.41	3.77	3.72	3.77	3.38	3.91	4.01	3.90	3.98	3.52	3.17	3.39	2.80	3.11	3.76	3.81	3.42	2.99	3.35	2.76];
  allN23    = [4.29	4.18	3.89	3.86	3.52	3.38	4.11	3.86	4.61	4.31	4.90	4.60	4.88	5.05	4.67	4.63	4.32	3.91	3.84	3.69	3.35	3.81	4.05	3.72	3.47	3.41	3.22];
  
%% Create TopoplotER  
  
load('/yourDirectory/Layouts/Layout_selected.mat');

data = struct;
data.freq = 1;
data.label = {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4',...
              'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','Cz','T8','P7',...
              'P3','Pz','P4','P8', 'O1','O2'}; % {1 x N} %spindle channels
data.dimord = 'chan_freq';
cfg = [];
cfg.layout = lay; %layy for spindle channels, lay for all channels with the rest set to 0
cfg.marker             = 'on';
cfg.colorbar           = 'yes';
cfg.interplimits       = 'head';%'head'; %'electrodes'; 
cfg.interpolation      = 'v4';
cfg.gridscale          = 250;
cfg.style              = 'straight'; %'both'; %'straight' 
cfg.colorbartext       = 'Spindle density (spindles / min)';
cfg.marker             = 'on'; %'on' %labels 
% cfg.markersize       = 3 %if marker = 'on'
cfg.markerfontsize     = 8; %8
cfg.markercolor        = [0 0 0];
cfg.markersize    = 24;
cfg.markersymbol    = '.';

%to higlight chosen channels
cfg.highlight          = 'marker'; %labels, marker
cfg.highlightchannel   =  {'C5','C3','C1','C2','C4','C6','CP3','CP4','FC3','FC4'};
cfg.highlightcolor = [1 1 1]; %uint8([240 33 61]);
cfg.highlightfontsize  = 10
cfg.highlightsize  = 24;
cfg.highlightsymbol  = '.'

cfg.xlim               = 'maxmin' ;
cfg.zlim               = 'maxmin'; %[4.6 18.1]
cfg.comment            = 'no';

% Cue
data.powspctrm = cueN23'; % [N x 1]
figure
cfg.zlim               = [3.09 5.34] %'maxmin'%[3.7 6.3];%[3.8 14] 
cue = ft_topoplotER(cfg,data)
c = sprintf('Cue period (N2 & N3)');
title({c,...
       '   '},'FontSize',16,'FontName','Calibri','FontWeight', 'normal');

   
% No cue
data.powspctrm = nocueN23'; % [N x 1]
figure
cfg.zlim               = [3.09 5.34]%'maxmin' % [10.29 18.1] 
nocue = ft_topoplotER(cfg,data)
nc = sprintf('No-cue period (N2 & N3)');
title({nc,...
       '   '},'FontSize',16,'FontName','Calibri', 'FontWeight', 'normal');

