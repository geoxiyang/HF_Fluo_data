% Plot SCOPE time-series run

work_dir = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/Kenia_2010_2015-04-05-1614/';

PSN     = dlmread([work_dir 'fluxes.dat'],'',[2,10,2920,10]);
t       = dlmread([work_dir 'fluxes.dat'],'',[2,3,2920,3]);
SIF760  = dlmread([work_dir 'fluorescence.dat'],'',[2,760-640+1,2920,760-640+1]);


% daily aggregation

doy     = 1:1:365;
SIF760_d= zeros(365,1);
PSN_d   = zeros(365,1);  

% QC
SIF760(SIF760<=0) = NaN;
PSN(PSN<=0)       = NaN;

for ii = 1:365
   
    SIF760_d(ii) = nanmean(SIF760(t>=ii & t<ii+1));
    PSN_d(ii)    = nanmean(PSN(t>=ii & t<ii+1));
    
end

%plotyy(1:1:365,SIF760_d,1:1:365,PSN_d); 


% monthly aggregation

month       = 1:1:12;
SIF760_m    = zeros(12,1);
PSN_m       = zeros(12,1);

for jj = 1:12
    
    st_m         = datenum(2010,jj,1) - datenum(2010,1,1) + 1;
    if jj == 12 
        en_m         = datenum(2011,1,1) - datenum(2010,1,1);
    else
        en_m         = datenum(2010,jj+1,1) - datenum(2010,1,1);
    end
    
    SIF760_m(jj) = nanmean(SIF760_d(doy>=st_m & doy<en_m));
    PSN_m(jj)    = nanmean(PSN_d(doy>=st_m & doy<en_m));
    
end

plotyy(1:1:12,SIF760_m,1:1:12,PSN_m);




