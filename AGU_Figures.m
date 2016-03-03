%% Figures for AGU
clear variables
clc


%% This is the figure for model vs. obs

work_dir1    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2013_2015-12-09-1326/';      % 2013 first half
work_dir2    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2014_2015-12-08-1152/';      % 2014 first half
work_dir3    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2014_2015-12-08-1745/';      % 2014 second half
work_dir4    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2013_2015-12-09-1749/';      % 2013 second half

data2013     = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013.mat','GPP','H','LE','doy_hourly','sif_hourly','doy_ems','apar_hourly','sunportion_hourly');
data2014     = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2014.mat','GPP','H','LE','doy_hourly','sif_hourly','doy_ems','apar_hourly','sunportion_hourly');


%% Read in data

fEnergy1     = dlmread([work_dir1 'fluorescence.dat'],'',2,0);
fEnergy2     = dlmread([work_dir2 'fluorescence.dat'],'',2,0);
fEnergy3     = dlmread([work_dir3 'fluorescence.dat'],'',2,0);
fEnergy4     = dlmread([work_dir4 'fluorescence.dat'],'',2,0);

fluxes1      = dlmread([work_dir1 'fluxes.dat'],'',2,0);
fluxes2      = dlmread([work_dir2 'fluxes.dat'],'',2,0);
fluxes3      = dlmread([work_dir3 'fluxes.dat'],'',2,0);
fluxes4      = dlmread([work_dir4 'fluxes.dat'],'',2,0);

% SCOPE simulations
LE_SCOPE1    = fluxes1(:,9);
LE_SCOPE2    = fluxes2(:,9);
LE_SCOPE3    = fluxes3(:,9);
LE_SCOPE4    = fluxes4(:,9);
A1           = fluxes1(:,11);
A2           = fluxes2(:,11);
A3           = fluxes3(:,11);
A4           = fluxes4(:,11);
f760_1       = fEnergy1(:,760-639);                     
f760_2       = fEnergy2(:,760-639);                     
f760_3       = fEnergy3(:,760-639);
f760_4       = fEnergy4(:,760-639);                     
time1        = fluxes1(:,4);
time2        = fluxes2(:,4);
time3        = fluxes3(:,4);
time4        = fluxes4(:,4);

aPAR2        = fluxes2(:,18);
aPAR3        = fluxes3(:,18);

A2           = fluxes2(:,11);
A3           = fluxes3(:,11);


doy_SCOPE1   = unique(fix(time1));
doy_SCOPE2   = unique(fix(time2))+365;
doy_SCOPE3   = unique(fix(time3))+365;
doy_SCOPE4   = unique(fix(time4));

doy_obs1     = unique(fix(data2013.doy_hourly));
doy_obs2     = unique(fix(data2014.doy_hourly))+365;

A1(A1<=0)     = NaN;
A2(A2<=0)     = NaN;
A3(A3<=0)     = NaN;
A4(A4<=0)     = NaN;
LE_SCOPE1(LE_SCOPE1<0) = NaN;
LE_SCOPE2(LE_SCOPE2<0) = NaN;
LE_SCOPE3(LE_SCOPE3<0) = NaN;
LE_SCOPE4(LE_SCOPE4<0) = NaN;
f760_1(f760_1<0.1) = NaN;
f760_2(f760_2<0.1) = NaN; 
f760_3(f760_3<0.1) = NaN; 
f760_4(f760_4<0.1) = NaN; 

%% Combine all the data

% A            = [A1;A2;A3];
% LE_SCOPE     = [LE_SCOPE1;LE_SCOPE2;LE_SCOPE3];
% f760         = [f760_1;f760_2;f760_3];
% doy_SCOPE    = [doy_SCOPE1;doy_SCOPE2;doy_SCOPE3];
% doy_obs      = [doy_obs1;doy_obs2];

f760_2014      = [f760_2;f760_3];
aPAR_2014      = [aPAR2;aPAR3];
A_2014         = [A2;A3];
time_2014      = [time2;time3];


plot(A_2014(time_2014>170 & time_2014<171)./aPAR_2014(time_2014>170 & time_2014<171),f760_2014(time_2014>170 & time_2014<171)./aPAR_2014(time_2014>170 & time_2014<171),'ro');


%% Plot


for ii = 1:length(doy_SCOPE1)
    
    lb = doy_SCOPE1(ii) + 0.25;
    ub = doy_SCOPE1(ii) + 0.75;
    daily_model_A1(ii)       = nanmean(A1(time1>=lb & time1<ub));
    daily_model_LE1(ii)      = nanmean(LE_SCOPE1(time1>=lb & time1<ub));
    daily_model_sif1(ii)     = nanmean(f760_1(time1>=lb & time1<ub));
 
end

for ii = 1:length(doy_SCOPE2)
    
    lb = doy_SCOPE2(ii) + 0.25;
    ub = doy_SCOPE2(ii) + 0.75;
    daily_model_A2(ii)       = nanmean(A2(time2+365>=lb & time2+365<ub));
    daily_model_LE2(ii)      = nanmean(LE_SCOPE2(time2+365>=lb & time2+365<ub));
    daily_model_sif2(ii)     = nanmean(f760_2(time2+365>=lb & time2+365<ub));
 
end

for ii = 1:length(doy_SCOPE3)
    
    lb = doy_SCOPE3(ii) + 0.25;
    ub = doy_SCOPE3(ii) + 0.75;
    daily_model_A3(ii)       = nanmean(A3(time3+365>=lb & time3+365<ub));
    daily_model_LE3(ii)      = nanmean(LE_SCOPE3(time3+365>=lb & time3+365<ub));
    daily_model_sif3(ii)     = nanmean(f760_3(time3+365>=lb & time3+365<ub));
 
end

for ii = 1:length(doy_SCOPE4)
    
    lb = doy_SCOPE4(ii) + 0.25;
    ub = doy_SCOPE4(ii) + 0.75;
    daily_model_A4(ii)       = nanmean(A4(time4>=lb & time4<ub));
    daily_model_LE4(ii)      = nanmean(LE_SCOPE4(time4>=lb & time4<ub));
    daily_model_sif4(ii)     = nanmean(f760_4(time4>=lb & time4<ub));
 
end

for jj = 1:length(doy_obs1)
   
    lb = doy_obs1(jj) + 0.25;
    ub = doy_obs1(jj) + 0.75;
    daily_obs_GPP1(jj) = nanmean(data2013.GPP(data2013.doy_hourly>=lb & data2013.doy_hourly<ub));
    daily_obs_LE1(jj)  = nanmean(data2013.LE(data2013.doy_hourly>=lb & data2013.doy_hourly<ub));
    daily_obs_sif1(jj) = nanmean(data2013.sif_hourly(data2013.doy_hourly>=lb & data2013.doy_hourly<ub));

end

for jj = 1:length(doy_obs2)
   
    lb = doy_obs2(jj) + 0.25;
    ub = doy_obs2(jj) + 0.75;
    daily_obs_GPP2(jj) = nanmean(data2014.GPP(data2014.doy_hourly+365>=lb & data2014.doy_hourly+365<ub));
    daily_obs_LE2(jj)  = nanmean(data2014.LE(data2014.doy_hourly+365>=lb & data2014.doy_hourly+365<ub));
    daily_obs_sif2(jj) = nanmean(data2014.sif_hourly(data2014.doy_hourly+365>=lb & data2014.doy_hourly+365<ub));

end

doy_SCOPE           = [doy_SCOPE1;doy_SCOPE4;doy_SCOPE2;doy_SCOPE3];
daily_model_A       = [daily_model_A1';daily_model_A4';daily_model_A2';daily_model_A3'];
daily_model_sif     = [daily_model_sif1';daily_model_sif4';daily_model_sif2';daily_model_sif3'];

doy_obs             = [doy_obs1;doy_obs2];
daily_obs_GPP       = [daily_obs_GPP1';daily_obs_GPP2'];
daily_obs_sif       = [daily_obs_sif1';daily_obs_sif2'];
sub_sif = daily_obs_sif >0 & daily_obs_sif <1;

fontsize = 16;
fontsize_small = 12;
fontname = 'Whitney';

subplot(2,1,1)
plot(doy_SCOPE,daily_model_A,'ro',doy_obs,daily_obs_GPP,'bo')
set(gca,'FontSize',fontsize,'FontName',fontname)
title('GPP','FontSize',fontsize,'FontName',fontname)
xlabel('Day of Year since 2013-01-01','FontSize',fontsize,'FontName',fontname)
ylabel({'GPP';'(umol/m^{2}/sec)'},'FontSize',fontsize,'FontName',fontname)
subplot(2,1,2)
plot(doy_SCOPE,daily_model_sif,'ro',doy_obs(sub_sif),daily_obs_sif(sub_sif),'bo')
set(gca,'FontSize',fontsize,'FontName',fontname)
title('SIF','FontSize',fontsize,'FontName',fontname)
xlabel('Day of Year since 2013-01-01','FontSize',fontsize,'FontName',fontname)
ylabel({'SIF';'(mw/m^{2}/nm/sr)'},'FontSize',fontsize,'FontName',fontname)

print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/model_obs_daily_comparison.png')
close(gcf)

%% The figure comparing GOSAT, GOME-2 and ground

% For GOME-2 pick 09:30 clear day; For GOSAT pick 13:30 clear day

% load('/Volumes/XiYangResearch/Projects/4.DiurnalLUE/2.Matlab/SIF.mat') 
% coordinate = [42.5,-72.2];
% 
% site_lat_s   = knnsearch(lat,coordinate(1));
% site_lon_s   = knnsearch(lon,coordinate(2));
% 
% site_lat     = lat(site_lat_s);
% site_lon     = lon(site_lon_s);
% 
% SIF_GOME_HF = SIF740(:,site_lon_s,site_lat_s).*0.582;
% for kk = 1:105
%     GOME_t(kk) = datenum(timeym(kk,1),timeym(kk,2),1) - datenum(2007,1,1);
% end

%plot(GOME_t,SIF_GOME_HF,'ro')
% 
% SIF_GOSAT_site_d            = zeros(62,1);
% SIF_GOSAT_1sigma_site_d     = zeros(62,1);
% SIF_GOME2_site_d            = zeros(62,1);
% SIF_GOME2_1sigma_site_d     = zeros(62,1);
% 
% for ii = 1:1 
%     
%     site_lat     = knnsearch(lat_gosat,coordinate(ii,1));
%     site_lon     = knnsearch(lon_gosat,coordinate(ii,2));
%     
%     site_lat_d_de(ii)= lat_gosat(site_lat);
%     site_lon_d_de(ii)= lon_gosat(site_lon);
%     
%     SIF_GOSAT_site_d(:,ii)        = SIF_GOSAT(site_lon,site_lat,:);
%     SIF_GOSAT_1sigma_site_d(:,ii) = SIF_1s_GOSAT(site_lon,site_lat,:);
%     SIF_GOME2_site_d(:,ii)        = SIF_GOME_760_RG(site_lon,site_lat,:);
%     SIF_GOME2_1sigma_site_d(:,ii) = SIF_GOME_760_1s_RG(site_lon,site_lat,:);
%     
% end




% Calculate monthly clear-day 09:30, 13:30 SIF


% data2013 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013.mat',...
%              'GPP',...
%              'H',...
%              'LE',...
%              'sif_hourly',...
%              'apar_hourly',...
%              'pri1_hourly',...
%              'pri2_hourly',...
%              'doy_hourly',...
%              'sunportion_hourly',...
%              'vpd_hourly',...
%              'ndvi_hourly',...
%              'evi_hourly');
%          
% data2014 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2014.mat',...
%              'GPP',...
%              'H',...
%              'LE',...
%              'sif_hourly',...
%              'apar_hourly',...
%              'pri1_hourly',...
%              'pri2_hourly',...
%              'doy_hourly',...
%              'sunportion_hourly',...
%              'vpd_hourly',...
%              'ndvi_hourly',...
%              'evi_hourly',...
%              'rin_hourly',...
%              'airT_hourly');
%          
%          
% temp_hour_2013      = data2013.doy_hourly - fix(data2013.doy_hourly);
% % gome_time_sub_2013  = find(temp_hour_2013 <= 9.6/24.0 & temp_hour_2013 >= 9.4/24.0);
% % gosat_time_sub_2013 = find(temp_hour_2013 <= 13.6/24.0 & temp_hour_2013 >= 13.4/24.0);
% % sp_flag_2013        = data2013.sunportion_hourly(gome_time_sub_2013) > 0.5;
% doy_2013            = fix(data2013.doy_hourly);
% 
% 
% temp_hour_2014      = data2014.doy_hourly - fix(data2014.doy_hourly);
% % gome_time_sub_2014  = find(temp_hour_2014 <= 9.6/24.0 & temp_hour_2014 >= 9.4/24.0);
% % gosat_time_sub_2014 = find(temp_hour_2014 <= 13.6/24.0 & temp_hour_2014 >= 13.4/24.0);
% % sp_flag_2014        = data2014.sunportion_hourly(gome_time_sub_2014) > 0.5;
% doy_2014            = fix(data2014.doy_hourly);
% 
% 
% for ii = 1:12
%     
%     lb = datenum(2013,ii,1) -  datenum(2013,1,1) + 1;
%     if ii == 12
%         ub = 365;
%     else
%         ub = datenum(2013,ii+1,1) -  datenum(2013,1,1);
%     end
%     
%     tower_gome_sif_2013(ii,1)   = nanmean(data2013.sif_hourly(temp_hour_2013 <= 9.6/24.0 & temp_hour_2013 >= 9.4/24.0 ...
%                                                         & data2013.doy_hourly >=lb & data2013.doy_hourly <=ub ...
%                                                         & data2013.sunportion_hourly >0.5));
%     tower_gosat_sif_2013(ii,1)  = nanmean(data2013.sif_hourly(temp_hour_2013 <= 13.6/24.0 & temp_hour_2013 >= 13.4/24.0 ...
%                                                         & data2013.doy_hourly >=lb & data2013.doy_hourly <=ub ...
%                                                         & data2013.sunportion_hourly >0.5));
%     date_gome_2013(ii,1)        = lb+datenum(2013,1,1)-datenum(2007,1,1);
%    
% end
% 
% for ii = 1:12
%     
%     lb = datenum(2014,ii,1) -  datenum(2014,1,1) + 1;
%     if ii == 12
%         ub = 365;
%     else
%         ub = datenum(2014,ii+1,1) -  datenum(2014,1,1);
%     end
%     
%     tower_gome_sif_2014(ii,1)   = nanmean(data2014.sif_hourly(temp_hour_2014 <= 10.1/24.0 & temp_hour_2014 >= 9.9/24.0 ...
%                                                         & data2014.doy_hourly >=lb & data2014.doy_hourly <=ub ...
%                                                         & data2014.sunportion_hourly >0.5));
%     tower_gosat_sif_2014(ii,1)  = nanmean(data2014.sif_hourly(temp_hour_2014 <= 13.6/24.0 & temp_hour_2014 >= 13.4/24.0 ...
%                                                         & data2014.doy_hourly >=lb & data2014.doy_hourly <=ub ...
%                                                         & data2014.sunportion_hourly >0.5));
%     date_gome_2014(ii,1)        = lb+datenum(2014,1,1)-datenum(2007,1,1);
%    
% end
% 
% tower_gome_sif = [tower_gome_sif_2013;tower_gome_sif_2014];
% doy_gome       = [date_gome_2013;date_gome_2014];
% 
% plot(GOME_t,SIF_GOME_HF,'ro',doy_gome,tower_gome_sif,'bo')
% 
% xlim([datenum(2013,6,1) - datenum(2007,1,1), datenum(2014,11,1) - datenum(2007,1,1)])







