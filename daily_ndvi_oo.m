% Calculate Daily mean NDVI from OceanOptic reflectance spectra

clear all
clc

load('SIF760_result.mat','raw_ref_result','wl');
load('hf_barn_2013_env.mat','cloud_ratio_daily');
ref = raw_ref_result(:,2:end);
[xdim_ref,ydim_ref] = size(ref);

ndvi_day      = zeros(130,3);
ndvi_day(:,1) = 170:1:299;

NIR_wl        = 750.0;
R_wl          = 685.0;  %705.0
NIR_ind       = wl >= NIR_wl & wl < NIR_wl+1;
R_ind         = wl >= R_wl   & wl < R_wl+1;

for ii = 1:xdim_ref
    NIR_temp = ref(ii,NIR_ind);
    R_temp  = ref(ii,R_ind);
    % Calculate NDVI at each sampling interval
    NIR(ii) = mean(NIR_temp(NIR_temp>=0.0 & NIR_temp<=1.0));
    R(ii)   = mean(R_temp(R_temp>=0.0 & R_temp<=1.0));
    ndvi(ii)= (NIR(ii) - R(ii))./(NIR(ii) + R(ii));
end

for uni_i = 6:6 %1:130
   
   % 1.Limit the data to the day of interest
   lb       = uni_i-1.0+170.;
   ub       = uni_i+170;
   sub_temp = raw_ref_result(:,1) >= lb & raw_ref_result(:,1) <= ub;
   
plot(raw_ref_result(sub_temp,1),ndvi(sub_temp),'ko')   
   
   % 2.Calculate Daily NDVI
   ndvi_tmp = ndvi(sub_temp);
   ndvi_day(uni_i,2) = mean(ndvi_tmp(ndvi_tmp>0));
   ndvi_day(uni_i,3) = std(ndvi_tmp(ndvi_tmp>0));     %STANDARD DEVIATION
    
end

plot(ndvi_day(:,1),ndvi_day(:,2),'ko');
ylim([0 1]);