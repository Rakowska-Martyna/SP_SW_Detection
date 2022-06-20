# SP_SW_Detection
MATLAB pipeline for **sleep spindles** and **slow waves detection** used in Rakowska et al. (2021)[^1] and (2022)[^2].

Contains 9 scripts for spindle analysis, numbered in the order in which they need to be run:
- **pr0** - create EEG layour for your data
- **pr1** - preprocesses raw EEG data for spindle analysis
- **pr2** - identifies spindle and slow-oscillation events across the sleep stage of interest & calculates relevant statistics (e.g. spindle density) 
- **pr3** - summarises spindle density from pr2 across participants 
- **pr4** - identifies spindles and slow-oscillation events without arousals
- **pr5** - computes the duration of the cue and no-cue period 
- **pr6**- identifies spindle and slow-oscillation events during the cue and no cue period of the chosen sleep stage. Calculates relevant statistics (e.g. spindle density).
- **pr7** - summarises spindle density from pr6 (during cue/no-cue period) across participants 
- **pr8** - creates topoplot figures (copied and pasted from Excel spreadsheet)

## Requirements 

The pipeline requires FieldTrip (https://www.fieldtriptoolbox.org/) and access to *myToolbox* by Miguel Navarrete (https://github.com/mnavarretem/myToolbox.git)

## Authors and contributors

* Martyna Rakowska
* Miguel Navarrete 

## Relevant methods section from Rakowska et al. (2021)[^1] :

*The relationship between spindles and behavioural outcomes was determined by focusing the analysis on 8 electrodes located over motor regions: 4 left (FC3, C5, C3, C1, CP3) and 4 right (FC4, C6, C4, C2, CP4). However, for visualisation purpose, the rest of the electrodes in the International 10–20 EEG system were also pre-processed and analysed as outlined below, and included in the final figure (Fig. 5A). Briefly, the raw data were first down-sampled to 250 Hz (for them to be comparable between the two EEG data acquisition systems) and filtered using a Chebyshev Type II infinite impulse response (IIR) filter (passband: f = [0.3 – 35] Hz; stopband: f 〈 0.1 Hz & f 〉 45 Hz). Then, for each participant, the channels were visually inspected and, if deemed noisy for the majority of the night, interpolated based on their triangulation-based neighbours. The final pre-processing step involved re-referencing the data to the mastoids (TP9, TP10). Algorithms for spindles and SOs counting (Navarrete et al., 2020) were subsequently employed to detect slow oscillations (0.3 – 2 Hz) and sleep spindles (11 – 16 Hz) at each electrode and in each sleep stage separately (N2, N3) or combined (N2 and N3). Briefly, for spindles detection, the data were filtered in a sigma band using an IIR filter again (passband: f = [11 – 16] Hz; stopband: f 〈 9 Hz & f 〉 18 Hz). Then, we used a 300 ms time window to compute the root mean squared (RMS) of the signal. Any event that had surpassed the 86.64 percentile (1.5 SD, Gaussian distribution) of the RMS signal was regarded as a candidate spindle. To fit the final spindle detection criteria (based on Iber et al. 2007), an event was deemed a sleep spindle if it occurred in the target sleep stage, lasted between 0.5 and 2.0 s and had at least 5 oscillations during that period (Navarrete et al., 2020). For SOs detection, the EEG data were filtered in the 0.3 – 2 Hz band using the IIR filter (passband: f = [0.3 – 2] Hz; stopband: f 〈 0.1 Hz & f 〉 4 Hz). Waves with negative deflection between −35 and −300 mV and with zero crossing between 0.13 and 1.66 s were considered SOs.* 

*The identified spindles and SOs were then separated into those that fell within the cue and no-cue periods. The cue period was defined as the time interval between 0 and 3.5 s after a tone onset (the longest inter-trial interval allowed), thus essentially encompassing the period from the onset of the first tone in a sequence until 3.5 s after the onset of the last one. The no-cue period was defined as the time interval between the sequences - from 3.5 to 20.0 s after the onset of the last tone in the sequence.*

*Spindle density was calculated by dividing the total number of spindles at each electrode by the length (in minutes) of the target period (cue period during target sleep stage, no-cue period during target sleep stage). Spindle density, together with the number of spindle and SO events during the cue and no-cue period of each of the target sleep stages, are presented in Table S2. Spindle laterality was obtained by subtracting spindle density over the right motor channels from the spindle density over the left motor channels.*

### References:

[^1]: Rakowska, M., Abdellahi, M. E., Bagrowska, P., Navarrete, M., & Lewis, P. A. (2021). Long term effects of cueing procedural memory reactivation during NREM sleep. NeuroImage, 244, 118573.
[^2]: Rakowska, M., Bagrowska, P., Lazari, A., Navarrete, M., Abdellahi, M. E., Johansen-Berg, H., & Lewis, P. A. (2022). Cueing motor memory reactivation during NREM sleep engenders learning-related changes in precuneus and sensorimotor structures. bioRxiv.
