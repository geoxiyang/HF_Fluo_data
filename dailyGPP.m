%% Calculate Daily Integral of GPP

doy = doy(doy>1);
gpp_raw = gpp_raw(doy>1);

gpp_doy = unique(floor(doy));
doy_size = size(gpp_doy);
gpp_day = zeros(doy_size(1),1);

% This calculate daily GPP
for ii = 1:doy_size(1)
   lb = double(gpp_doy(ii))+0.;
   ub = double(gpp_doy(ii))+1.;
   gpp_temp = gpp_raw(doy>= lb & doy <ub);
   %Convert umol m-2 s-1 to g C m-2 day-1
   gpp_day(ii) = sum(gpp_temp*12*1e-6*30*60);       
end

% This calculate 8-day GPP corresponds to MODIS GPP
% From DOY 169-297
% Calculate accumulated GPP every 8 days, for example 169 - 177

accu_gpp = zeros(16,1);
for jj = 1:16
    accu_gpp(jj) = sum(gpp_day(gpp_doy>= 169 + (jj-1)*8 & gpp_doy< 169 + jj*8));   
end

save('HF_2013_GPP.mat');
