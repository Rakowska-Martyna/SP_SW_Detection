# SP_SW_Detection
MATLAB pipeline for sleep spindles and slow waves detection

Contains 12 scripts for spindle analysis, numbered in the order in which they need to be run:
- pr1 - preprocesses raw EEG data for spindle analysis
- pr2 - identifies spindle and slow-oscillation events across the sleep stage of interest & calculates relevant statistics (e.g. spindle density) 
- pr3 - summarises spindle density from pr2 across participants 
- pr4 - identifies spindles and slow-oscillation events without arousals
- pr5 - computes the duration of the cue and no-cue period 
- pr6- identifies spindle and slow-oscillation events during the cue and no cue period of the chosen sleep stage. Calculates relevant statistics (e.g. spindle density).
- pr7 - summarises spindle density from pr6 (during cue/no-cue period) across participants 
- pr8 - creates topoplot figures (copied and pasted from Excel spreadsheet)

The pipeline requires access to myToolbox by Miguel Navarrete available at https://github.com/mnavarretem/myToolbox.git
