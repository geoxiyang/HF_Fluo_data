% Compare SIF calculated using SFM and SVD

%% Section 1. Compare diurnal patterns
load('hf_2013_svd.mat','final_svd','final_sfm');
% load('SIF760_result.mat','final_sfm');
load('hf_2013_svd.mat','rad_time','irrad_time');

day = 175;

svd_time_idx = find(final_svd(:,1) >= day+0.2 & final_svd(:,1) <= day+0.8 & final_svd(:,2) < 3.0);
sfm_time_idx = find(final_sfm(:,1) >= day+0.2 & final_sfm(:,1) <= day+0.8 & final_sfm(:,3) > 0.9);

subplot(3,1,1)
plot(final_svd(svd_time_idx,1),final_svd(svd_time_idx,2),'ro',final_sfm(sfm_time_idx,1),final_sfm(sfm_time_idx,2),'bo');
subplot(3,1,2)
plot(final_svd(svd_time_idx,1),final_svd(svd_time_idx,3),'ko');
subplot(3,1,3)
irrad_time_1 = sortrows(irrad_time,1);
plot(final_svd(svd_time_idx,1),1,'ro')

%% Section 2. Compare seasonal pattern
% 
% svd = load('SIF760daily_2013_SVD.mat','SIF_mean');
% sfm = load('SIF760daily.mat','SIF_mean');
% 
% plot(svd.SIF_mean(:,1),svd.SIF_mean(:,2),'ro',sfm.SIF_mean(:,1),sfm.SIF_mean(:,2),'bo')