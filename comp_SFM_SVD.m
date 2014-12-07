% Compare SIF calculated using SFM and SVD

%% Section 1. Compare diurnal patterns
load('hf_2013_svd.mat','final_svd');
load('SIF760_result.mat','final_result_time');

% day = 199;
% 
% svd_time_idx = find(final_svd(:,1) >= day+0.2 & final_svd(:,1) <= day+0.8);
% sfm_time_idx = find(final_result_time(:,1) >= day+0.2 & final_result_time(:,1) <= day+0.8);
% 
% subplot(2,1,1)
% plot(final_svd(svd_time_idx,1),final_svd(svd_time_idx,2),'ro',final_result_time(sfm_time_idx,1),final_result_time(sfm_time_idx,2),'bo');
% subplot(2,1,2)
% plot(final_svd(svd_time_idx,1),final_svd(svd_time_idx,3),'ko');

%% Section 2. Compare seasonal pattern

svd = load('SIF760daily_2013_SVD.mat','SIF_mean');
sfm = load('SIF760daily.mat','SIF_mean');

plot(svd.SIF_mean(:,1),svd.SIF_mean(:,2),'ro',sfm.SIF_mean(:,1),sfm.SIF_mean(:,2),'bo')