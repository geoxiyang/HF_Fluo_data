%% Check individual days data, visualization etc.

clear variables
clc

datapath = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';
data2012 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2012.mat');      
data2013 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newcutoff.mat'); 
data2014 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2014.mat');
% raw2012  = load([datapath 'SIF760_result_2012_newtest.mat']);
% raw2013  = load([datapath 'SIF760_result_2013_newtest.mat']);
% raw2014  = load([datapath 'SIF760_result_2014.mat']);

% Convention: 2012: Red; 2013: Blue; 2014: Black
% doy_tmp_2012    = unique(fix(data2012.doy_hourly));
% doy_tmp_2013    = unique(fix(data2013.doy_hourly));
% doy_tmp_2014    = unique(fix(data2014.doy_hourly));
% 
% figure
% plot(doy_tmp_2012,data2012.gpp_daily_2012,'ro')
% hold on
% plot(doy_tmp_2013,data2013.gpp_daily_2013,'bo')
% plot(doy_tmp_2014,data2014.gpp_daily_2014,'ko')
% hold off

% figure
% plot(raw2012.final_result_time(raw2012.final_result_time(:,3)>0.9,1),raw2012.final_result_time(raw2012.final_result_time(:,3)>0.9,2)*0.6,'ro')
% hold on
% plot(raw2013.final_result_time(raw2013.final_result_time(:,3)>0.9,1),raw2013.final_result_time(raw2013.final_result_time(:,3)>0.9,2),'bo')
% plot(raw2014.final_result_time(raw2014.final_result_time(:,3)>0.9,1),raw2014.final_result_time(raw2014.final_result_time(:,3)>0.9,2),'ko')
% hold off
% ylim([0 5])

% figure
% plot(data2012.doy_ems,data2012.GPP,'ro')
% hold on
% plot(data2013.doy_ems,data2013.GPP,'bo')
% plot(data2014.doy_ems,data2014.GPP,'ko')
% hold off

figure
plot(doy_tmp_2012,sif_daily_2012,'ro')
hold on
plot(doy_tmp_2013,sif_daily_2013,'bo')
plot(doy_tmp_2014,sif_daily_2014,'ko')
hold off