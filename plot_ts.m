% Plot SCOPE time-series run

work_dir = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2014-09-17-1338/';

PSN     = dlmread([work_dir 'fluxes.dat'],'',[2,10,4121,10]);
t       = dlmread([work_dir 'fluxes.dat'],'',[2,3,4121,3]);
SIF760  = dlmread([work_dir 'fluorescence.dat'],'',[2,760-640+1,4121,760-640+1]);

