%============
%  SIFyield vs VPD analysis
%  Xi Yang, geoxiyang@gmail.com
%============

clear
clc

datapath = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';
load([datapath,'SIF760daily.mat'],'halfhourly_result');
load([datapath,'HF_2013_GPP.mat']);
load([datapath,'hf_barn_2013_env.mat'],'apar','VPD','cloud_ratio');

sif_cube = reshape(halfhourly_result(:,2),48,130);
gpp_cube = reshape(gpp_raw,48,130);
apar_cube= reshape(apar,48,130);
vpd_cube = reshape(VPD,48,130);
cr_cube  = reshape(cloud_ratio,48,130);

sif_cube(sif_cube<=0 | sif_cube>=3) = NaN;
gpp_cube(gpp_cube<=0) = NaN;
apar_cube(apar_cube<=0) = NaN;
vpd_cube(vpd_cube<=0) = NaN;
cr_cube(cr_cube>1) = NaN;

% calculate SIF yield and LUE

sify_cube = sif_cube./apar_cube;
lue_cube  = gpp_cube./apar_cube;

ave_indice  = [13,13.5,14]*2;
sunny_sub   = bsxfun(@gt,nansum(cr_cube(ave_indice,:) < 0.5,1),0) ;


sify_mean   = nanmean(sify_cube(ave_indice,sunny_sub),1);
vpd_mean    = nanmean(vpd_cube(ave_indice,sunny_sub),1);
lue_mean    = nanmean(lue_cube(ave_indice,:),1);
sif_mean    = nanmean(sif_cube(ave_indice,:),1);
apar_mean   = nanmean(apar_cube(ave_indice,:),1);
gpp_mean    = nanmean(gpp_cube(ave_indice,:),1);
plot(vpd_mean,sify_mean,'mo','MarkerSize',12)



% ylim([0 0.0015])

% plot vpd vs. sify and lue at a specific time of day throughout the season

% time_indice = 13.5*2;
% 
% plot(sify_cube(time_indice,:),vpd_cube(time_indice,:),'ro','MarkerSize',12)