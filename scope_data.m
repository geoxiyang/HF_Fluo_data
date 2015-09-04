%% Make data files for SCOPE runs
% Since HF Fluxes (final results) are in hourly format, we will need to
% recalculate SIF. In the end, we will need a table with the following
% columns (output filenames in parenthesis). We will need to separate the
% data for year 2013 and 2014.
% 1. time (t_.dat) in decimal day of year
% 2. year (year_.dat)
% 3. TOC incoming shortwave rad (Rin_.dat) w m-2
% 4. TOC incoming longwave rad  (Rli_.dat) w m-2
% 5. Air pressure (p_.dat) hPa or mbar
% 6. Air temperature above the canopy (Ta_.dat) Celsius
% 7. Vapor pressure above the canopy (ea_.dat) hPa or mbar
% 8. Windspeed (u_.dat) m s-1
%
% A separate file for hourly:
% 1. GPP (umol m-2 s-1)
% 2. SIF (mw m-2 sr-1 s-1)
% 3. LE  (w m-2)
% 4. H   (w m-2)
%
% In addition, we calculate Vapor Pressure Deficit (VPD):
% VPD = ea * (100-RH)/RH
% ========================================================================
% Note:
% TOC incoming shortwave, longwave radiation, air temperature are from 
% Andrew's Barn tower dataset. Time in the set is for the end of the half hour.
% Air pressure and windspeed are from Bill's EMS tower data, time is for
% the beginning for the hour.
% GPP, LE, H are from Bill's EMS tower data.
% SIF is from Xi's FluoSpec, recalculated in this code to get hourly data.
%
% Xi Yang (geoxiyang@gmail.com). Aug.4, 2015.

clc
clear variables
%% Section 1: Calculate hourly data

%  read data
radfilepath         = '/Volumes/XiYangResearch/Data/6.HFData/hf249-01-radiometric1.xlsx';
datapath            = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';
fluxesfilepath1     = '/Volumes/XiYangResearch/Data/6.HFData/hf004-02-filled.xlsx';
%fluxesfilepath1     = '/Volumes/XiYangResearch/Data/6.HFData/US-Ha1-2014-Results.xlsx';
fluxesfilepath2     = '/Volumes/XiYangResearch/Data/6.HFData/hf004-01-final.xlsx';
%fluxesfilepath2     = '/Volumes/XiYangResearch/Data/6.HFData/EC1.14.xlsx'; %2014 windspeed
airPfilepath        = '/Volumes/XiYangResearch/Data/6.HFData/hf001-10-15min-m.xlsx'; %There are some missing air pressure data from EMS, so we have to use this one

% for 2012
% load([datapath,'SIF760daily_2012_newtest.mat'],'halfhourly_result');
% [rawdata,txt,raw]   = xlsread(radfilepath,'B12123:BK15914');
% [rawflux,txt2,raw2] = xlsread(fluxesfilepath1,'A182018:AG183913');
% [rawflux2,txt3,raw3]= xlsread(fluxesfilepath2,'A182018:AL183913');
% [rawflux3,txt4,raw4]= xlsread(airPfilepath,'Q265921:Q273600');

% % for 2013
load([datapath,'SIF760daily_2013_newcutoff.mat'],'halfhourly_result');
[rawdata,txt,raw]   = xlsread(radfilepath,'B27531:BK33770');
[rawflux,txt2,raw2] = xlsread(fluxesfilepath1,'A189698:AG192817');
[rawflux2,txt3,raw3]= xlsread(fluxesfilepath2,'A189698:AL192817');
[rawflux3,txt4,raw4]= xlsread(airPfilepath,'Q296737:Q309216');

% % for 2014
% load([datapath,'SIF760daily_2014_newtest.mat'],'halfhourly_result');
% [rawdata,txt,raw]   = xlsread(radfilepath,'B42987:BK51290');
% [rawflux,txt2,raw2] = xlsread(fluxesfilepath1,'A3027:BE7178');
% [rawflux2,txt3,raw3]= xlsread(fluxesfilepath2,'A3026:AL7177');
% [rawflux3,txt4,raw4]= xlsread(airPfilepath,'Q327650:Q344257');

doy                 = rawdata(:,1);
rin_30min           = rawdata(:,4);
rli_30min           = rawdata(:,10);
% 2012 Andrew's airt and rh data are not there (NA, comment out)
airT_30min          = rawdata(:,61);
rh_30min            = rawdata(:,62);
tot_ppfd_30min      = rawdata(:,16);
ref_ppfd_30min      = rawdata(:,17);
blw1_ppfd_30min     = rawdata(:,18);
blw2_ppfd_30min     = rawdata(:,19);
blw3_ppfd_30min     = rawdata(:,20);
blw4_ppfd_30min     = rawdata(:,21);
apar_ppfd_30min     = tot_ppfd_30min - ref_ppfd_30min - nanmean([blw1_ppfd_30min,blw2_ppfd_30min,blw3_ppfd_30min,blw4_ppfd_30min],2);
pri1_30min          = rawdata(:,53);
pri2_30min          = rawdata(:,55);
direct_rad          = rawdata(:,24);
total_rad           = rawdata(:,22);
sunny_portion       = direct_rad./total_rad;

A                   = 8.07131;  
B                   = 1730.63;
C                   = 233.426;
% For 2013 and 2014
ea_30min            = 10.^(A-B./(C+airT_30min));
vpd_30min           = ea_30min .* (100-rh_30min)./rh_30min;

GPP                 = abs(rawflux(:,13));
%GPP                 = rawflux(:,47);  %2014
H                   = rawflux2(:,16);
%H                   = rawflux(:,6);   %2014
%airT_ems            = rawflux(:,8);   %2014
airT_ems            = rawflux(:,17);    %rawflux2(:,28);
rh_ems              = rawflux2(:,23);
LE                  = rawflux2(:,17).*((2.502*1e6-(2.308*10e3.*airT_ems))*18*1e-3); %*1e-3; % NOTE: 2012 LE needs to multiply 1e-3 probably because the original data's unit is mol/m2
% LE                  = rawflux(:,49); % 2014 used original
doy_ems             = rawflux2(:,4);
% doy_ems             = rawflux2(:,1); %2014
windspeed           = rawflux2(:,5);
% windspeed           = rawflux2(:,2); %2014

airP                = rawflux3(:,1);
% 2012
% doy_airP            = (215+1/96):1/96:294;
% 2013
 doy_airP            = (170+1/96):1/96:300;
% 2014
% doy_airP            = 127:1/96:(300-1/96);

% calculate hourly data
% 2013: 130; 2014: 173; 2012: 79
totdays             = 130;
doy_hourly          = nan(totdays*24,1);
rin_hourly          = nan(totdays*24,1);
rli_hourly          = nan(totdays*24,1);
airT_hourly         = nan(totdays*24,1);
rh_hourly           = nan(totdays*24,1);
ea_hourly           = nan(totdays*24,1);
sif_hourly          = nan(totdays*24,1);
airP_hourly         = nan(totdays*24,1);
apar_hourly         = nan(totdays*24,1);
pri1_hourly         = nan(totdays*24,1);
pri2_hourly         = nan(totdays*24,1);
vpd_hourly          = nan(totdays*24,1);
sunportion_hourly   = nan(totdays*24,1);
% for 2012 airt and rh are hourly
% airT_hourly         = airT_ems;
% rh_hourly           = rh_ems;
% ea_hourly           = 10.^(A-B./(C+airT_ems));

GPP(GPP<=0) = NaN;
LE(LE<=0)   = NaN;

for ii = 1:totdays
   for jj = 1:24
       
         lb = double(ii+doy(1)-1)+(jj-1)/24;
         ub = double(ii+doy(1)-1)+jj/24;
         
         doy_hourly((ii-1)*24+jj)   = lb;
         rin_hourly((ii-1)*24+jj)   = nanmean(rin_30min(doy>lb & doy<=ub));
         rli_hourly((ii-1)*24+jj)   = nanmean(rli_30min(doy>lb & doy<=ub));
         apar_hourly((ii-1)*24+jj)  = nanmean(apar_ppfd_30min(doy>lb & doy<=ub));
         pri1_hourly((ii-1)*24+jj)  = nanmean(pri1_30min(doy>lb & doy<=ub));
         pri2_hourly((ii-1)*24+jj)  = nanmean(pri2_30min(doy>lb & doy<=ub));
         
        airT_hourly((ii-1)*24+jj)  = nanmean(airT_30min(doy>lb & doy<=ub));
        rh_hourly((ii-1)*24+jj)    = nanmean(rh_30min(doy>lb & doy<=ub));
        ea_hourly((ii-1)*24+jj)    = nanmean(ea_30min(doy>lb & doy<=ub));
        vpd_hourly((ii-1)*24+jj)   = nanmean(vpd_30min(doy>lb & doy<=ub));
         
         sif_hourly((ii-1)*24+jj)   = nanmean(halfhourly_result(halfhourly_result(:,1)>=lb & halfhourly_result(:,1)<ub,2));
         
         airP_hourly((ii-1)*24+jj)  = nanmean(airP(doy_airP>lb & doy_airP<=ub));
         
         sunportion_hourly((ii-1)*24+jj) = nanmean(sunny_portion(doy>lb & doy<=ub));
   end
end

vpd_hourly           = ea_hourly .* (100-rh_hourly)./rh_hourly;
save('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newcutoff.mat');




%% Section 2: Arrange dataset for outputs

% load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly.mat');
% output_folder = '/Volumes/XiYangResearch/src/SCOPE/data/input/dataset HF_ts/';
% 
% year = repmat(2013,[length(doy_hourly),1]);
% rin_hourly(rin_hourly<0) = 0.0;
% 
% dlmwrite([output_folder 'year_.dat'],year);
% dlmwrite([output_folder 't_.dat'],doy_hourly);
% dlmwrite([output_folder 'Rin_.dat'],rin_hourly);
% dlmwrite([output_folder 'Rli_.dat'],rli_hourly);
% dlmwrite([output_folder 'p_.dat'],airP_hourly);
% dlmwrite([output_folder 'Ta_.dat'],airT_hourly);
% dlmwrite([output_folder 'ea_.dat'],ea_hourly);
% dlmwrite([output_folder 'u_.dat'],windspeed);






























