%% Calculate Daily Integral of GPP, LE, and H

% doy = doy(doy>1);
% gpp_raw = gpp_raw(doy>1);
% 
% gpp_doy = unique(floor(doy));
% doy_size = size(gpp_doy);
% gpp_day = zeros(doy_size(1),1);
% 
% % This calculate daily GPP
% for ii = 1:doy_size(1)
%    lb = double(gpp_doy(ii))+0.;
%    ub = double(gpp_doy(ii))+1.;
%    gpp_temp = gpp_raw(doy>= lb & doy <ub);
%    %Convert umol m-2 s-1 to g C m-2 day-1
%    gpp_day(ii) = sum(gpp_temp*12*1e-6*30*60);       
% end
% 
% % This calculate 8-day GPP corresponds to MODIS GPP
% % From DOY 169-297
% % Calculate accumulated GPP every 8 days, for example 169 - 177
% 
% accu_gpp = zeros(16,1);
% for jj = 1:16
%     accu_gpp(jj) = sum(gpp_day(gpp_doy>= 169 + (jj-1)*8 & gpp_doy< 169 + jj*8));   
% end
% 
% save('HF_2013_GPP.mat');


%% For 2014
% load ?US-Ha1-2014-Results.txt?
uniq_doy = unique(DoY);
doy_size = size(uniq_doy);

%gpp_day  = zeros(doy_size(1),1);
le_day   = zeros(doy_size(1),1);
%gpp_day_m= zeros(doy_size(1),1);
le_day_m = zeros(doy_size(1),1);
%h_day    = zeros(doy_size(1),1);

%GPP_f(GPP_f<=0) = nan;
LE_orig(LE_orig<=0)   = nan;

for ii = 1:doy_size(1)
   lb = double(uniq_doy(ii))+0.;
   ub = double(uniq_doy(ii))+1.;
   
%   gpp_temp = GPP_f(DoY>= lb & DoY <ub);
   le_temp  = LE_orig(DoY>= lb & DoY <ub);
   %Convert umol m-2 s-1 to g C m-2 day-1 -- harvard forest 2014 data is
   %hourly
%   gpp_day(ii)      = nansum(gpp_temp*12*1e-6*60*60);
%   gpp_day_m(ii)    = nanmean(gpp_temp*12*1e-6*60*60);
   le_day(ii)       = nansum(le_temp);
   le_day_m(ii)     = nanmean(le_temp);
     
end

