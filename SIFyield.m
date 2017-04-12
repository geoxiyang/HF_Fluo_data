%% === New SIFyield analysis ======================%

data2012 = load('/Volumes/XiYangBackUp/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2012.mat',...
             'GPP',...
             'H',...
             'LE',...
             'sif_hourly',...
             'apar_hourly',...
             'pri1_hourly',...
             'pri2_hourly',...
             'doy_hourly',...
             'sunportion_hourly',...
             'vpd_hourly',...
             'ndvi_hourly',...
             'evi_hourly');
         
data2013 = load('/Volumes/XiYangBackUp/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013.mat',...
             'GPP',...
             'H',...
             'LE',...
             'sif_hourly',...
             'apar_hourly',...
             'pri1_hourly',...
             'pri2_hourly',...
             'doy_hourly',...
             'sunportion_hourly',...
             'vpd_hourly',...
             'ndvi_hourly',...
             'evi_hourly');
         
data2014 = load('/Volumes/XiYangBackUp/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2014_newtimestamp.mat',...
             'GPP',...
             'H',...
             'LE',...
             'sif_hourly',...
             'apar_hourly',...
             'pri1_hourly',...
             'pri2_hourly',...
             'doy_hourly',...
             'sunportion_hourly',...
             'vpd_hourly',...
             'ndvi_hourly',...
             'evi_hourly',...
             'rin_hourly',...
             'airT_hourly');

data2012.sif_hourly(data2012.sif_hourly<=0) = NaN;
data2013.sif_hourly(data2013.sif_hourly<=0) = NaN;
data2014.sif_hourly(data2014.sif_hourly<=0) = NaN;         

% ====== Calculate SZA ==============
for kk = 1: length(data2014.doy_hourly)
    
    dnumber      = data2014.doy_hourly(kk) + datenum(2014,1,1);
    dvec         = datevec(dnumber);
    time1.year   = 2014;
    time1.month  = dvec(2);
    time1.day    = dvec(3);
    time1.hour   = dvec(4);
    time1.min    = dvec(5);
    time1.sec    = dvec(6);
    time1.UTC    = -5;
    location.latitude = 42.5;
    location.longitude = -72.2;
    location.altitude = 100; 
    
    sun_pos = sun_position(time1,location);
    sza(kk)  = sun_pos.zenith;
    sun_azimuth(kk) = sun_pos.azimuth;
end



hours     = data2014.doy_hourly - fix(data2014.doy_hourly);



th        = 0.5;       
sify_2014 = data2014.sif_hourly./data2014.apar_hourly;
lue_2014  = data2014.GPP./data2014.apar_hourly;

scatter(lue_2014,sify_2014,'ro')
xlim([0,0.05])
ylim([0,0.0011])
sub1 = data2014.sunportion_hourly > 0.6;% & hours > 0.582 & hours < 0.584;% & data2014.apar_hourly > 1000; % & data2014.doy_hourly >170 & data2014.doy_hourly <171;% & data2014.apar_hourly < 600 & data2014.apar_hourly > 100;% & data2014.vpd_hourly <30;% & sza'< 30;% & sza'< 60;
%scatter(lue_2014(sub1),sify_2014(sub1),[], hours(sub1),'filled')

scatter(sza(sub1),data2014.sif_hourly(sub1),[],hours(sub1),'filled')

%scatter(lue_2014(sub1),sify_2014(sub1),[],sza(sub1),'filled'); %data2014.evi_hourly(sub1)
% scatter(data2014.GPP(sub1),data2014.sif_hourly(sub1),[],sza(sub1),'filled')
% hold on
% sub2 = data2014.sunportion_hourly > 0.6 & data2014.apar_hourly > 0 & sza'< 60 & sza'> 30;
% scatter(lue_2014(sub2),sify_2014(sub2),[],sza(sub2),'filled')
% hold on
% sub3 = data2014.sunportion_hourly > 0.6 & data2014.apar_hourly > 0 & sza'> 60;
% scatter(lue_2014(sub3),sify_2014(sub3),[],sza(sub3),'filled')
%  
% hold off
% lsline
colorbar
 
 
% colorbar
% figure
% plot(lue_2014(sub),sify_2014(sub),'ro')
% figure
% scatter(lue_2014(sub),sify_2014(sub),[],data2014.evi_hourly(sub),'filled')
% colorbar
% scatter(lue_2014(sub),sify_2014(sub),[],data2014.apar_hourly(sub),'filled')
% colorbar

% subplot(2,1,1)
% scatter(data2014.apar_hourly,lue_2014,[],data2014.sunportion_hourly,'filled')
% colorbar
% ylim([0 0.1])
% xlabel('APAR')
% ylabel('LUE')
% 
% subplot(2,1,2)
% plot(data2014.apar_hourly(sub),lue_2014(sub),'ro')

% print('-dpng','-r300','/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/GPP_SIF_ColorbySZA_ThresholdSunportion06_2014.png')  %APAR0SZALARGER60
% close(gcf)












