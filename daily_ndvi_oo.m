% Calculate Daily mean NDVI from OceanOptic reflectance spectra



% clear variables
% clc

load('SIF760_result_2014.mat','raw_ref_result','wl');
%load('hf_barn_2013_env.mat','cloud_ratio_daily');
ref = raw_ref_result(:,2:end);
[xdim_ref,ydim_ref] = size(ref);
time = raw_ref_result(:,1);

ndvi_day      = zeros(173,3); %130,3
ndvi_day(:,1) = 127:1:299;
chl_day(:,1)  = 127:1:299;
NIR_wl        = 750.0;
R_wl          = 705.0;
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
     ndvi(ndvi<=0) = NaN;
%     chl = 40.65*ndvi.^2+121.88.*ndvi-0.77;

for uni_i = 1:173
   
   % 1.Limit the data to the day of interest
   lb       = uni_i-1.0+127;
   ub       = uni_i+127;
   sub_temp = raw_ref_result(:,1) >= lb & raw_ref_result(:,1) <= ub;
   
%    plot(raw_ref_result(sub_temp,1)-fix(raw_ref_result(sub_temp,1)),ndvi(sub_temp),'ko')   
%    ylim([0.8 0.9])
%    hold on
   
   
%    % 2.Calculate Daily NDVI
   ndvi_tmp = ndvi(sub_temp);
   %chl_tmp  = chl(sub_temp);
   ndvi_day(uni_i,2) = mean(ndvi_tmp(ndvi_tmp>0));
   ndvi_day(uni_i,3) = std(ndvi_tmp(ndvi_tmp>0));     %STANDARD DEVIATION
end

    chl_day(:,2) = ((30-10)/(max(ndvi_day(:,2))-min(ndvi_day(:,2)))).*(ndvi_day(:,2)-min(ndvi_day(:,2)))+10;     %40 and 10 are from my RSE manuscript %90 10 are from Dillen etal. % NOTE: 90 as the max Vcmax seems to be too high see Bonan JGR 2012

%hold off
plot(chl_day(:,1),chl_day(:,2),'ko');
%ylim([0 1]);