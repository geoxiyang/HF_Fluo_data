% interpolate 30 minutes Harvard Barn Tower data to 15 minutes
% Also we used these data to calculate daily values

% !!!!! Windows and Mac read xls differently

[data_30min,txt,raw] = xlsread('/Volumes/XiYangResearch/Projects/9.Fluorescence/7.WeatherData/Harvard Barn Radiometric Data Master 12 Feb 2014.xlsx' ...
                        ,'B27534:DB33773');  %'A27726:DB27774'
                    
load('HF_2013_GPP.mat','gpp_raw');
load('SIF760daily.mat','halfhourly_result')
doy = data_30min(:,1);
%Temp = data_30min(:,104);
Rshort = data_30min(:,5);
% ndvi = data_30min(:,94);
% 
% f655 = data_30min(:,85);
% f860 = data_30min(:,86);
% f470 = data_30min(:,87);
% f530 = data_30min(:,83);
% f570 = data_30min(:,84);
% out_570 = data_30min(:,76);
% in_570  = data_30min(:,68);
% out_530 = data_30min(:,75);
% in_530  = data_30min(:,67);
% 
% 
% f530_1 = out_530./in_530;
% f570_1 = out_570./in_570;
% 
% diffuse_rad = data_30min(:,43);
% direct_rad = data_30min(:,44);
% total_rad = data_30min(:,42);
% 
below1_rad = data_30min(:,41);
below2_rad = data_30min(:,40);
below3_rad = data_30min(:,19);
ref_rad = data_30min(:,18);
above_rad = data_30min(:,17);
% 
% cloud_ratio = diffuse_rad ./total_rad;
% apar = above_rad - ref_rad - nanmean([below1_rad,below2_rad,below3_rad],2);
% fapar = apar./above_rad;
% par = above_rad;
% evi = 2.5 * (f860-f655)./(f860 + 6*f655-7.5*f470+1);
% pri = (f530-f570)./(f530+f570);
% pri2= data_30min(:,91);
% pri3= (f530_1-f570_1)./(f530_1+f570_1);
% gpp_raw(gpp_raw<0) = nan;
% apar(apar<0) = nan;
% lue_raw = gpp_raw./apar;
% lue_yield = lue_raw./(halfhourly_result(:,2)./apar);


%UNIQUE DOY
unidoy = unique(floor(doy));
Rshort_daily =  nan(numel(unidoy),1);
% ndvi_daily = nan(numel(unidoy),1);
% evi_daily = nan(numel(unidoy),1);
% pri_daily = nan(numel(unidoy),1);
% pri2_daily= nan(numel(unidoy),1);
% pri3_daily= nan(numel(unidoy),1);
% % LAI = nan(numel(unidoy),1);
% cloud_ratio_daily = nan(numel(unidoy),1);
% apar_daily = nan(numel(unidoy),1);
% fapar_daily = nan(numel(unidoy),1);
% par_daily = nan(numel(unidoy),1);
% apar_sum = nan(numel(unidoy),1);
% lue_daily = nan(numel(unidoy),1);
% lue_midday = nan(numel(unidoy),1);
% lue_yield_midday = nan(numel(unidoy),1);
%tmin_daily = nan(numel(unidoy),1);
for kk = 1:numel(unidoy)  %1:numel(unidoy)
   % Calculate PRI daily value
   lb = double(unidoy(kk))+0/24;
   ub = double(unidoy(kk))+24/24;
   
   Rshort_daily(kk) = mean(Rshort(doy>=lb & doy<ub));
   
%   tmin_daily(kk) = min(Temp(doy>=lb & doy<=ub));
%    lb_mid = double(unidoy(kk))+11/24;
%    ub_mid = double(unidoy(kk))+13/24;
%    apar_daily(kk) = mean(apar(doy>=lb & doy<=ub));
%    fapar_daily(kk) = mean(fapar(doy>=lb & doy<=ub));
%    par_daily(kk) = mean(par(doy>=lb & doy<=ub));
%     lue_daily(kk) = nanmean(lue_raw(doy>=lb & doy<=ub));
%     lue_midday(kk) = nanmean(lue_raw(doy>=lb_mid & doy<=ub_mid));
%     lue_yield_midday(kk) = nanmean(lue_yield(doy>=lb_mid & doy<=ub_mid));
% %    apar_sum(kk) = sum(apar(doy>=lb & doy<=ub)*30*60);
%   cloud_ratio_daily(kk) = mean(cloud_ratio(doy>=lb & doy<=ub));
%    
%    % PRI is the average of that day
%     pri_daily(kk) = nanmean(pri(doy>=lb & doy<ub & abs(pri) <= 0.2));
%     pri2_daily(kk)= nanmean(pri2(doy>=lb & doy<ub & abs(pri) <= 0.2));
%     pri3_daily(kk)= nanmean(pri2(doy>=lb & doy<ub & abs(pri) <= 0.2));
%   ndvi_daily(kk) = nanmean(ndvi(doy>=lb & doy<ub));
%    evi_daily(kk) = nanmean(evi(doy>=lb & doy<ub));
%    
%    % LAI calculated in the following way:
%    % First, gap fraction P = below canopy PAR/above canopy PAR
%    % Below canopy PAR is the average of two below canopy measurements
%    below1_temp = below1_rad(doy>=lb & doy<ub);
%    below2_temp = below2_rad(doy>=lb & doy<ub);
%    total_temp = total_rad(doy>=lb & doy<ub);
%    gap_p = zeros(numel(total_temp),1);
%    gap_p(:,1) = (below1_temp(:,1) + below2_temp(:,1))./total_temp(:,1);
%    % Second, find the gap_p value with the solar zenith angle (SZA)=57
%    % Need to calculate the SZA for each time the measurement was made
%    doy_temp = doy(doy>=lb & doy<ub);
%    sza = zeros(numel(doy_temp),1);
%    for jj = 1:numel(doy_temp)
%      [year,month,day,hour,min,sec] = datevec(datenum('1/1/2013')+doy(jj));
%      time = struct('year',year,'month',month,'day',day,'hour',hour,...
%                    'min',min,'sec',sec,'UTC',-5);
%      location = struct('longitude',-72.1899,'latitude',42.5353,'altitude',350);
%      sun_loc = sun_position(time,location);
%      sza(jj,1) = sun_loc.zenith;
%    end
%    sza_morning_sub = knnsearch(sza(1:numel(sza)/2,1),57,'K',1);
%    sza_after_sub = knnsearch(sza((numel(sza)/2):end,1),57,'K',1);
%    % Calculate an AM and a PM LAI, and then average them
%    LAI1 = -log(gap_p(sza_morning_sub))/(0.5/cos(1)); % 1 radian = 57.2 degree
%    LAI2 = -log(gap_p(sza_after_sub+numel(sza)/2-1))/(0.5/cos(1));
%    
%    LAI(kk) = (LAI1 + LAI2)/2;
%    
   %
   
   
end



% par_daily(par_daily<0) = NaN;
% apar_daily(apar_daily<0) = NaN;
% fapar_daily(apar_daily<0) = NaN;
save('hf_barn_2013_env.mat','-append')

% total_rad = data_30min(:,42);
% rshort_dn = data_30min(:,6);
% rlong_dn = data_30min(:,12);

% t15 = doy(1)+(0:95)*15/(24*60);
% rshort_dn_15 = zeros(96,1);
% rlong_dn_15 = zeros(96,1);
% total_rad_hfems = NaN(6128,1);
% rh_hfems = NaN(6128,1);
% 
% 
% 
% 
% for ii = 1:6240
%     
%     
%     
%     
% %     if (mod(ii,2) == 0)
%       [subs,dists] = knnsearch(doy_hum(:,1),hfems_date(ii,1),'K',2);
%       rh_hfems(ii,1) = (rh(subs(1),1) + rh(subs(2),1))/2.0;
%   
% %       [subs,dists] = knnsearch(doy(:,1),t15(1,ii),'K',2);
% %       rshort_dn_15(ii,1) = (rshort_dn(subs(1),1) + rshort_dn(subs(2),1))/2.0;
% %       rlong_dn_15(ii,1) = (rlong_dn(subs(1),1) + rlong_dn(subs(2),1))/2.0;      
% %     else
% %       rshort_dn_15(ii,1) = rshort_dn(floor(ii/2)+1,1);
% %       rlong_dn_15(ii,1) = rlong_dn(floor(ii/2)+1,1);
% %     end
% end

