%% make figures for SIF_SCOPE_comparison

clear variable
clc

data2012 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2012.mat',...
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
         
data2013 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013.mat',...
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
         
data2014 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2014.mat',...
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
%          
% need to make daily SIF, GPP, LE, APAR, PRI. Here we use data of all kinds of
% sky conditions

% QC

data2012.sif_hourly(data2012.sif_hourly<=0) = NaN;
data2013.sif_hourly(data2013.sif_hourly<=0) = NaN;
data2014.sif_hourly(data2014.sif_hourly<=0) = NaN;


doy_tmp_2012    = unique(fix(data2012.doy_hourly));
doy_tmp_2013    = unique(fix(data2013.doy_hourly));
doy_tmp_2014    = unique(fix(data2014.doy_hourly));

gpp_daily_2012  = nan(length(doy_tmp_2012),1);
gpp_daily_2013  = nan(length(doy_tmp_2013),1);
gpp_daily_2014  = nan(length(doy_tmp_2014),1);
le_daily_2012   = nan(length(doy_tmp_2012),1);
le_daily_2013   = nan(length(doy_tmp_2013),1);
le_daily_2014   = nan(length(doy_tmp_2014),1);
apar_daily_2012 = nan(length(doy_tmp_2012),1);
apar_daily_2013 = nan(length(doy_tmp_2013),1);
apar_daily_2014 = nan(length(doy_tmp_2014),1);
sif_daily_2012  = nan(length(doy_tmp_2012),1);
sif_daily_2013  = nan(length(doy_tmp_2013),1);
sif_daily_2014  = nan(length(doy_tmp_2014),1);
% EF is only calculated for daytime (6am to 6pm)
ef_daily_2012   = nan(length(doy_tmp_2012),1);
ef_daily_2013   = nan(length(doy_tmp_2013),1);
ef_daily_2014   = nan(length(doy_tmp_2014),1); 

ndvi_daily_2012   = nan(length(doy_tmp_2012),1);
ndvi_daily_2013   = nan(length(doy_tmp_2013),1);
ndvi_daily_2014   = nan(length(doy_tmp_2014),1); 
evi_daily_2012   = nan(length(doy_tmp_2012),1);
evi_daily_2013   = nan(length(doy_tmp_2013),1);
evi_daily_2014   = nan(length(doy_tmp_2014),1);

rin_daily_2014      = nan(length(doy_tmp_2014),1);
airT_daily_2014     = nan(length(doy_tmp_2014),1);

gpp_th_daily_2013  = nan(length(doy_tmp_2013),1);
gpp_th_daily_2014  = nan(length(doy_tmp_2014),1);
apar_th_daily_2013 = nan(length(doy_tmp_2013),1);
apar_th_daily_2014 = nan(length(doy_tmp_2014),1);
sif_th_daily_2013  = nan(length(doy_tmp_2013),1);
sif_th_daily_2014  = nan(length(doy_tmp_2014),1);
threshold          = 0.5;

for ii = 1:length(doy_tmp_2012)
    
    gpp_daily_2012(ii,1)    = nanmean(data2012.GPP(data2012.GPP >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<(doy_tmp_2012(ii)+1)));
    le_daily_2012(ii,1)     = nanmean(data2012.LE(data2012.LE >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    sif_daily_2012(ii,1)    = nanmean(data2012.sif_hourly(data2012.sif_hourly >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    apar_daily_2012(ii,1)   = nanmean(data2012.apar_hourly(data2012.apar_hourly>0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    vpd_daily_2012(ii,1)    = nanmean(data2012.vpd_hourly(data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    ef_daily_2012(ii,1)     = nanmean(data2012.LE(data2012.LE >0 & data2012.H >0 & data2012.doy_hourly>=doy_tmp_2012(ii)+0.25 & data2012.doy_hourly<doy_tmp_2012(ii)+0.75)./...
                                     (data2012.H(data2012.LE >0 & data2012.H >0 & data2012.doy_hourly>=doy_tmp_2012(ii)+0.25 & data2012.doy_hourly<doy_tmp_2012(ii)+0.75)+...
                                     data2012.LE(data2012.LE >0 & data2012.H >0 & data2012.doy_hourly>=doy_tmp_2012(ii)+0.25 & data2012.doy_hourly<doy_tmp_2012(ii)+0.75)));
    ndvi_daily_2012(ii,1)   = nanmean(data2012.ndvi_hourly(data2012.ndvi_hourly >0 & data2012.doy_hourly>=(doy_tmp_2012(ii)+0.4) & data2012.doy_hourly<(doy_tmp_2012(ii)+0.6)));
    evi_daily_2012(ii,1)    = nanmean(data2012.evi_hourly(data2012.evi_hourly >0 & data2012.doy_hourly>=(doy_tmp_2012(ii)+0.4) & data2012.doy_hourly<(doy_tmp_2012(ii)+0.6)));
    
    
end

for ii = 1:length(doy_tmp_2013)
    
    gpp_daily_2013(ii,1)    = nanmean(data2013.GPP(data2013.GPP >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    le_daily_2013(ii,1)     = nanmean(data2013.LE(data2013.LE >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    sif_daily_2013(ii,1)    = nanmean(data2013.sif_hourly(data2013.sif_hourly >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    apar_daily_2013(ii,1)   = nanmean(data2013.apar_hourly(data2013.apar_hourly>0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    vpd_daily_2013(ii,1)    = nanmean(data2013.vpd_hourly(data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    ef_daily_2013(ii,1)     = nanmean(data2013.LE(data2013.LE >0 & data2013.H >0 & data2013.doy_hourly>=doy_tmp_2013(ii)+0.25 & data2013.doy_hourly<doy_tmp_2013(ii)+0.75)./...
                                     (data2013.H(data2013.LE >0 & data2013.H >0 & data2013.doy_hourly>=doy_tmp_2013(ii)+0.25 & data2013.doy_hourly<doy_tmp_2013(ii)+0.75)+...
                                     data2013.LE(data2013.LE >0 & data2013.H >0 & data2013.doy_hourly>=doy_tmp_2013(ii)+0.25 & data2013.doy_hourly<doy_tmp_2013(ii)+0.75)));
    ndvi_daily_2013(ii,1)   = nanmean(data2013.ndvi_hourly(data2013.ndvi_hourly >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    evi_daily_2013(ii,1)    = nanmean(data2013.evi_hourly(data2013.evi_hourly >0 & data2013.doy_hourly>=(doy_tmp_2013(ii)+0.4) & data2013.doy_hourly<(doy_tmp_2013(ii)+0.6)));
    
    gpp_th_daily_2013(ii,1) = nanmean(data2013.GPP(data2013.GPP >0 & data2013.sunportion_hourly > threshold & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    sif_th_daily_2013(ii,1) = nanmean(data2013.sif_hourly(data2013.sif_hourly >0 & data2013.sunportion_hourly > threshold & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    apar_th_daily_2013(ii,1)= nanmean(data2013.apar_hourly(data2013.apar_hourly>0 & data2013.sunportion_hourly > threshold & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    
    
end

for ii = 1:length(doy_tmp_2014)
    
    gpp_daily_2014(ii,1)    = nanmean(data2014.GPP(data2014.GPP >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    le_daily_2014(ii,1)     = nanmean(data2014.LE(data2014.LE >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    sif_daily_2014(ii,1)    = nanmean(data2014.sif_hourly(data2014.sif_hourly >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    apar_daily_2014(ii,1)   = nanmean(data2014.apar_hourly(data2014.apar_hourly>0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    vpd_daily_2014(ii,1)    = nanmean(data2014.vpd_hourly(data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    ef_daily_2014(ii,1)     = nanmean(data2014.LE(data2014.LE >0 & data2014.H >0 & data2014.doy_hourly>=doy_tmp_2014(ii)+0.25 & data2014.doy_hourly<doy_tmp_2014(ii)+0.75)./...
                                     (data2014.H(data2014.LE >0 & data2014.H >0 & data2014.doy_hourly>=doy_tmp_2014(ii)+0.25 & data2014.doy_hourly<doy_tmp_2014(ii)+0.75)+...
                                     data2014.LE(data2014.LE >0 & data2014.H >0 & data2014.doy_hourly>=doy_tmp_2014(ii)+0.25 & data2014.doy_hourly<doy_tmp_2014(ii)+0.75)));    
    ndvi_daily_2014(ii,1)   = nanmean(data2014.ndvi_hourly(data2014.ndvi_hourly >0 & data2014.doy_hourly>=(doy_tmp_2014(ii)+0.4) & data2014.doy_hourly<(doy_tmp_2014(ii)+0.6)));
    evi_daily_2014(ii,1)    = nanmean(data2014.evi_hourly(data2014.evi_hourly >0 & data2014.doy_hourly>=(doy_tmp_2014(ii)+0.4) & data2014.doy_hourly<(doy_tmp_2014(ii)+0.6)));      

    airT_daily_2014(ii,1)   = nanmean(data2014.airT_hourly(data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    rin_daily_2014(ii,1)    = nanmean(data2014.rin_hourly(data2014.rin_hourly >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    
    gpp_th_daily_2014(ii,1) = nanmean(data2014.GPP(data2014.GPP >0 & data2014.sunportion_hourly > threshold &data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    sif_th_daily_2014(ii,1) = nanmean(data2014.sif_hourly(data2014.sif_hourly >0 & data2014.sunportion_hourly > threshold & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    apar_th_daily_2014(ii,1)= nanmean(data2014.apar_hourly(data2014.apar_hourly>0 & data2014.sunportion_hourly > threshold & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    
end

gpp_daily_all   = [gpp_daily_2012;gpp_daily_2013;gpp_daily_2014];
le_daily_all    = [le_daily_2012;le_daily_2013;le_daily_2014];
apar_daily_all  = [apar_daily_2012;apar_daily_2013;apar_daily_2014];
sif_daily_all   = [sif_daily_2012;sif_daily_2013;sif_daily_2014];
doy_daily_all   = [doy_tmp_2012;doy_tmp_2013+366;doy_tmp_2014+366+365];
vpd_daily_all   = [vpd_daily_2012;vpd_daily_2013;vpd_daily_2014];
ef_daily_all    = [ef_daily_2012;ef_daily_2013;ef_daily_2014];
ndvi_daily_all  = [ndvi_daily_2012;ndvi_daily_2013;ndvi_daily_2014];
evi_daily_all   = [evi_daily_2012;evi_daily_2013;evi_daily_2014];

gpp_th_daily_all = [gpp_th_daily_2013;gpp_th_daily_2014];
apar_th_daily_all= [apar_th_daily_2013;apar_th_daily_2014];
sif_th_daily_all = [sif_th_daily_2013;sif_th_daily_2014];


gpp_hourly_all   = [data2012.GPP;data2013.GPP;data2014.GPP];
le_hourly_all    = [data2012.LE;data2013.LE;data2014.LE];
apar_hourly_all  = [data2012.apar_hourly;data2013.apar_hourly;data2014.apar_hourly];
sif_hourly_all   = [data2012.sif_hourly;data2013.sif_hourly;data2014.sif_hourly];
doy_hourly_all   = [data2012.doy_hourly;data2013.doy_hourly+366;data2014.doy_hourly+366+365];
sp_hourly_all    = [data2012.sunportion_hourly;data2013.sunportion_hourly;data2014.sunportion_hourly];

gpp_th_hourly_all= [data2013.GPP;data2014.GPP];
apar_th_hourly_all=[data2013.apar_hourly;data2014.apar_hourly];
sif_th_hourly_all=[data2013.sif_hourly;data2014.sif_hourly];
sp_th_hourly_all= [data2013.vpd_hourly;data2014.vpd_hourly];
vpd_th_hourly_all=[data2013.vpd_hourly;data2014.vpd_hourly];

gpp_th_hourly_all(gpp_th_hourly_all<0|sp_th_hourly_all<threshold) = NaN;
apar_th_hourly_all(apar_th_hourly_all<0|sp_th_hourly_all<threshold) = NaN;
sif_th_hourly_all(sif_th_hourly_all<0|sp_th_hourly_all<threshold) = NaN;

sif_daily_all(sif_daily_all>1.1) = NaN;


% subplot(2,1,1)
% plot([doy_tmp_2013+366;doy_tmp_2014+366+365],gpp_th_daily_all,'ro')
% subplot(2,1,2)
% plot([doy_tmp_2013+366;doy_tmp_2014+366+365],sif_th_daily_all,'ko')

% th = 0.5;
% 
% gpp_hourly_all(gpp_hourly_all<=0 | sp_hourly_all<th)     = NaN;
% le_hourly_all(le_hourly_all<=0 | sp_hourly_all<th)       = NaN;
% apar_hourly_all(apar_hourly_all<=0 | sp_hourly_all<th)   = NaN;
% sif_hourly_all(sif_hourly_all<=0 | sp_hourly_all<th)     = NaN;


%% SZA fixed

% % calculate SZA for each hour, using 2013 as a start
% 
% 
% for kk = 1: length(data2013.doy_hourly)
%     
%     dnumber      = data2013.doy_hourly(kk) + datenum(2013,1,1);
%     dvec         = datevec(dnumber);
%     time1.year   = 2013;
%     time1.month  = dvec(2);
%     time1.day    = dvec(3);
%     time1.hour   = dvec(4);
%     time1.min    = dvec(5);
%     time1.sec    = dvec(6);
%     time1.UTC    = -5;
%     location.latitude = 42.5;
%     location.longitude = -72.2;
%     location.altitude = 100; 
%     
%     sun_pos = sun_position(time1,location);
%     sza(kk)  = sun_pos.zenith;
% 
% 
% end
% 
% 
% sub = sza>=0 & sza<=25;
% 
% plot(data2013.sif_hourly(sub),data2013.GPP(sub),'ro')


%% New code making Fig.1

% sif_daily_all(sif_daily_all>1.1) = NaN;
% 
% fontsize = 16;
% fonttype = 'Gill Sans MT';
% 
% %set(gcf,'Color',[0 0 0]);
% 
% subplot(4,1,2)
% plot(doy_daily_all,gpp_daily_all,'co');
% hold on
% line([365,365],[0,40],'Color','w','LineStyle','--','LineWidth',2)
% line([731,731],[0,40],'Color','w','LineStyle','--','LineWidth',2)
% text(50,25,num2str(2012),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(415,25,num2str(2013),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(780,25,num2str(2014),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% hold off
% set(gca,'Color',[0 0 0],...
%     'FontSize',fontsize,...
%     'FontName',fonttype,...
%     'XLim',[0,1100],...
%     'YLim',[0,30],...
%     'XColor',[1 1 1],...
%     'YColor',[1 1 1])
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel({'GPP';'(umol/m^{2}/sec)'},'FontSize',fontsize,'FontName',fonttype);
% 
% % subplot(5,1,5)
% % plot(doy_daily_all,le_daily_all,'bo')
% % set(gca,'FontSize',fontsize,'YLim',[0, 200],'FontName',fonttype)
% % %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% % ylabel('LE(w/m^{2})','FontSize',fontsize,'FontName',fonttype);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% subplot(4,1,3)
% plot(doy_daily_all,evi_daily_all,'g.','MarkerSize',30)
% hold on
% line([365,365],[0.3,0.8],'Color','w','LineStyle','--','LineWidth',2)
% line([731,731],[0.3,0.8],'Color','w','LineStyle','--','LineWidth',2)
% text(50,(25/30)*0.8,num2str(2012),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(415,(25/30)*0.8,num2str(2013),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(780,(25/30)*0.8,num2str(2014),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% 
% hold off
% set(gca,'Color',[0 0 0],...
%     'FontSize',fontsize,...
%     'YLim',[0.3,0.8],...
%     'XLim',[0,1100],...
%     'FontName',fonttype,...
%     'XColor',[1 1 1],...
%     'YColor',[1 1 1])
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel('EVI','FontSize',fontsize,'FontName',fonttype);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% subplot(4,1,1)
% plot(doy_daily_all,sif_daily_all,'ro')
% hold on
% line([365,365],[0.0,1.2],'Color','w','LineStyle','--','LineWidth',2)
% line([731,731],[0.0,1.2],'Color','w','LineStyle','--','LineWidth',2)
% text(50,(25/30)*1.2,num2str(2012),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(415,(25/30)*1.2,num2str(2013),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(780,(25/30)*1.2,num2str(2014),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% 
% hold off
% set(gca,'Color',[0 0 0],...
%     'FontSize',fontsize,...
%     'FontName',fonttype,...
%     'YLim',[0.0,1.2],...
%     'XLim',[0,1100],...
%     'XColor',[1 1 1],...
%     'YColor',[1 1 1])
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel({'SIF';'(mw/m^{2}/nm/sr)'},'FontSize',fontsize);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% subplot(4,1,4)
% plot(doy_daily_all,apar_daily_all,'o','Color',[1 .5 0])
% hold on
% line([365,365],[0,900],'Color','k','LineStyle','--','LineWidth',2,'Color','w')
% line([731,731],[0,900],'Color','k','LineStyle','--','LineWidth',2,'Color','w')
% text(50,(25/30)*900,num2str(2012),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(415,(25/30)*900,num2str(2013),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% text(780,(25/30)*900,num2str(2014),'FontSize',fontsize,'FontName',fonttype,'Color','w')
% 
% hold off
% set(gca,'Color',[0 0 0],...
%     'FontSize',fontsize,...
%     'FontName',fonttype,...
%     'YLim',[0.0,900],...
%     'XLim',[0,1100],...
%     'XColor',[1 1 1],...
%     'YColor',[1 1 1])
% xlabel('Day of Year since 2012-01-01','FontSize',fontsize,'FontName',fonttype);
% ylabel({'APAR';'(umol/m^{2}/sec)'},'FontSize',fontsize,'FontName',fonttype);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% export_fig('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIF_GPP_APAR_HF.png','-r300','-dpng','-transparent','-nocrop')

%print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIF_GPP_APAR_HF.png')
%close(gcf)


%% Fig.2 

% %---------------Fig.2a--------------------
% %-----------------------------------------
% figure
% 
% 
% 

% FontName = 'Whitney';
% 
% plot(sif_daily_all,gpp_daily_all,'r.','MarkerSize',20);
% %Label axes
% % LE(w/m^{2})%GPP(umol/m^{2}/sec)
% 
% hold on
% % Michaelis-Menten fit
% ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Robust = 'Off';
% opts.StartPoint = [max(gpp_daily_all) 0.5];
% opts.Lower  = [0 0];
% opts.Upper  = [2*max(gpp_daily_all) 5];
% 
% % Fit model to data.
% [fitresult, gof] = fit( sif_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all)), gpp_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all)), ft, opts );
% 
% % h = plot(fitresult);
% % set(h,'Color','k','LineWidth',2)
% % set(gca,'YLim',[0,max(gpp_daily_all)*1.2],...
% %     'FontSize',16);
% % legend('off')
% % 
% % hh = text(0.05,max(gpp_daily_all)*1.1,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')],'FontName',FontName);
% % hh.FontSize = 20;   
% 
% % Linear fit
% ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts2.Display = 'Off';
% opts2.Robust = 'Off';
% opts2.StartPoint = [20 1];
% opts2.Lower  = [0 0];
% opts2.Upper  = [2*max(gpp_daily_all) 5];
% 
% % Fit model to data.
% [fitresult2, gof2] = fit( sif_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all)), gpp_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all)), ft2, opts2);
% 
% h = plot(fitresult2);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(gpp_daily_all)*1.2],...
%     'FontSize',16,'FontName',FontName);
% legend('off')
% 
% hh = text(0.05,max(gpp_daily_all)*1.0,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')],'FontName',FontName);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20,'FontName',FontName);
% ylabel('GPP(umol/m^{2}/sec)','FontSize',20,'FontName',FontName);
% hold off
% 
% set(gca,'FontName',FontName)
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_SIF_GPP_linearonly.png')
% close(gcf)
% % %-----------------------------------------
% % %---------------Fig.2b--------------------
% % %-----------------------------------------
% figure
% plot(sif_daily_all,le_daily_all,'b.','MarkerSize',20);
% % Label axes
% %LE(w/m^{2})%GPP(umol/m^{2}/sec)
% 
% hold on
% % Michaelis-Menten fit
% ft = fittype('a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Robust = 'Off';
% opts.StartPoint = [max(le_daily_all) 0.5];
% opts.Lower  = [0 0];
% opts.Upper  = [2*max(le_daily_all) 5];
% 
% % Fit model to data.
% [fitresult, gof] = fit( sif_daily_all(~isnan(sif_daily_all) & ~isnan(le_daily_all)), le_daily_all(~isnan(sif_daily_all) & ~isnan(le_daily_all)), ft, opts );
% 
% h = plot(fitresult);
% set(h,'Color','k','LineWidth',2)
% set(gca,'YLim',[0,max(le_daily_all)*1.05],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,185,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
% hh.FontSize = 20;   
% 
% % Linear fit
% ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts2.Display = 'Off';
% opts2.Robust = 'Off';
% opts2.StartPoint = [20 1];
% opts2.Lower  = [0 0];
% opts2.Upper  = [2*max(le_daily_all) 5];
% 
% % Fit model to data.
% [fitresult2, gof2] = fit( sif_daily_all(~isnan(sif_daily_all) & ~isnan(le_daily_all)), le_daily_all(~isnan(sif_daily_all) & ~isnan(le_daily_all)), ft2, opts2);
% 
% h = plot(fitresult2);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(le_daily_all)*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,170,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('LE(w/m^{2})','FontSize',20);
% 
% set(gca,'YLim',[0,200],...
%     'FontSize',16);
% hold off
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_SIF_LE.png')
% close(gcf)
% 
% %-----------------------------------------
% %---------------Fig.2c--------------------
% %-----------------------------------------
% figure
% plot(apar_daily_all(~isnan(sif_daily_all)& ~isnan(gpp_daily_all) & ~isnan(apar_daily_all)),...
%     sif_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all) & ~isnan(apar_daily_all)),...
%     'm.','MarkerSize',20);
% hold on
% 
% % Linear fit
% ft3 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts3.Display = 'Off';
% opts3.Robust = 'Off';
% opts3.StartPoint = [0.001 0.005];
% opts3.Lower  = [-Inf -Inf];
% opts3.Upper  = [Inf Inf];
% 
% testa = sif_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all) & ~isnan(apar_daily_all)); %sif_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all));
% testb = apar_daily_all(~isnan(sif_daily_all) & ~isnan(gpp_daily_all) & ~isnan(apar_daily_all) );
% % Fit model to data.
% [fitresult3, gof3] = fit(testb, testa, ft3, opts3);
% 
% h = plot(fitresult3);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'XLim',[0,max(testb)*1.2],...
%     'YLim',[0 1],...
%     'FontSize',16);
% legend('off')
% 
% %0.05
% hh = text(10,0.9,['Linear fit r^{2}=' num2str(gof3.rsquare,'%4.2f')],'FontName',FontName);
% hh.FontSize = 20; 
% 
% ylabel('SIF(mw/m^{2}/nm/sr)','FontSize',20,'FontName',FontName); %'GPP(umol/m^{2}/sec)',
% xlabel('APAR(umol/m^{2}/sec)','FontSize',20,'FontName',FontName); %
% hold off
% 
% set(gca,'FontName',FontName)
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_dailymean_APAR_SIF.png')
% close(gcf)

%% Fig.2d-f

% 


% filename1        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_GPP.png';
% filename2        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_LE.png';
% filename3        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_APAR.png';
% 
% %---------------Fig.2d--------------------
% %-----------------------------------------
% figure
% 
% % temp_hours = [data2012.doy_hourly-fix(data2012.doy_hourly);data2013.doy_hourly-fix(data2013.doy_hourly);data2014.doy_hourly-fix(data2014.doy_hourly)];
% % temp_subs  = temp_hours>=9.9/24.0 & temp_hours<=10.1/24.0;
% % 
% % sif_hourly_all = sif_hourly_all(temp_subs);
% % gpp_hourly_all = gpp_hourly_all(temp_subs);
% % 
% plot(sif_hourly_all,gpp_hourly_all,'r.','MarkerSize',20);
% 
% 
% 
% 
% % Label axes
% %LE(w/m^{2})%GPP(umol/m^{2}/sec)
% 
% hold on
% % Michaelis-Menten fit
% ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Robust = 'Off';
% opts.StartPoint = [max(gpp_hourly_all) 0.5];
% opts.Lower  = [0 0];
% opts.Upper  = [2*max(gpp_hourly_all) 5];
% 
% % Fit model to data.
% [fitresult, gof] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), gpp_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), ft, opts );
% 
% h = plot(fitresult);
% set(h,'Color','r','LineWidth',2)
% set(gca,'YLim',[0,max(gpp_hourly_all)*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,max(gpp_hourly_all)*1.1,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
% hh.FontSize = 20;   
% 
% % Linear fit
% ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts2.Display = 'Off';
% opts2.Robust = 'Off';
% opts2.StartPoint = [20 1];
% opts2.Lower  = [0 0];
% opts2.Upper  = [2*max(gpp_hourly_all) 5];
% 
% % Fit model to data.
% [fitresult2, gof2] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), gpp_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), ft2, opts2);
% 
% h = plot(fitresult2);
% set(h,'Color','r','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(gpp_hourly_all)*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,max(gpp_hourly_all)*1.0,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('GPP(umol/m^{2}/sec)','FontSize',20);
% hold off
% 
% print('-dpng','-r300',filename1)
% close(gcf)
% %-----------------------------------------
% %---------------Fig.2b--------------------
% %-----------------------------------------
% figure
% plot(sif_hourly_all,le_hourly_all,'b.','MarkerSize',20);
% % Label axes
% %LE(w/m^{2})%GPP(umol/m^{2}/sec)
% 
% hold on
% % Michaelis-Menten fit
% ft = fittype('a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Robust = 'Off';
% opts.StartPoint = [max(le_hourly_all) 0.5];
% opts.Lower  = [0 0];
% opts.Upper  = [2*max(le_hourly_all) 5];
% 
% % Fit model to data.
% [fitresult, gof] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), le_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), ft, opts );
% 
% h = plot(fitresult);
% set(h,'Color','k','LineWidth',2)
% set(gca,'YLim',[0,max(le_hourly_all)*1.05],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,500,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
% hh.FontSize = 20;   
% 
% % Linear fit
% ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts2.Display = 'Off';
% opts2.Robust = 'Off';
% opts2.StartPoint = [20 1];
% opts2.Lower  = [0 0];
% opts2.Upper  = [2*max(le_hourly_all) 5];
% 
% % Fit model to data.
% [fitresult2, gof2] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), le_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), ft2, opts2);
% 
% h = plot(fitresult2);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(le_hourly_all)*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,460,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('LE(w/m^{2})','FontSize',20);
% 
% set(gca,'YLim',[0,max(le_hourly_all)*1.05],...
%     'FontSize',16);
% hold off
% 
% print('-dpng','-r300',filename2)
% close(gcf)
% 
% %-----------------------------------------
% %---------------Fig.2c--------------------
% %-----------------------------------------
% figure
% plot(sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)),apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)),'m.','MarkerSize',20);
% 
% hold on
% 
% % Linear fit
% ft3 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts3.Display = 'Off';
% opts3.Robust = 'Off';
% opts3.StartPoint = [500 120];
% opts3.Lower  = [-Inf -Inf];
% opts3.Upper  = [Inf Inf];
% 
% testa = sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all));
% testb = apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all));
% % Fit model to data.
% [fitresult3, gof3] = fit(testa, testb, ft3, opts3);
% 
% h = plot(fitresult3);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)))*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,2200,['Linear fit r^{2}=' num2str(gof3.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('APAR(umol/m^{2}/sec)','FontSize',20);
% hold off
% 
% print('-dpng','-r300',filename3)
% close(gcf)

%% Fig.3 Monthly diurnal comparison

% for month_i = 1:36
%     
%    lb = datenum(2012,month_i,1) - datenum(2012,1,1) +1;
%    ub = datenum(2012,month_i+8,1) - datenum(2012,1,1);
% 
%    temp_sif   = sif_hourly(doy_hourly>=lb & doy_hourly<ub);
%    temp_apar  = apar_hourly(doy_hourly>=lb & doy_hourly<ub);
%    temp_GPP   = GPP(doy_hourly>=lb & doy_hourly<ub);
%    temp_LE    = LE(doy_hourly>=lb & doy_hourly<ub);
%    temp_sp    = sunportion_hourly(doy_hourly>=lb & doy_hourly<ub);
%    temp_doy   = doy_hourly(doy_hourly>=lb & doy_hourly<ub);
%    
%    temp_doy1  = unique(floor(temp_doy));
% 
% 
% end

%% Fig.4 
% 
% lue_daily_all   = gpp_daily_all./apar_daily_all;
% sify_daily_all  = sif_daily_all./apar_daily_all;
% lue_hourly_all  = gpp_hourly_all./apar_hourly_all;
% sify_hourly_all = sif_hourly_all./apar_hourly_all;
% 
% lue_th_daily_all = gpp_th_daily_all./apar_th_daily_all;
% sify_th_daily_all= sif_th_daily_all./apar_th_daily_all;
% 
% lue_th_hourly_all = gpp_th_hourly_all./apar_th_hourly_all;
% sify_th_hourly_all= sif_th_hourly_all./apar_th_hourly_all;
% 
% % scatter(lue_th_hourly_all,sify_th_hourly_all,24,vpd_th_hourly_all,'filled')
% 
% 
% 
% sub              = lue_daily_all<0.5;
% sub1             = lue_daily_all<0.5 & apar_daily_all>500;
% sub2             = lue_daily_all<0.5 & apar_daily_all<500;
% lue_daily_all1   = lue_daily_all(sub1);
% sify_daily_all1  = sify_daily_all(sub1);
% lue_daily_all2   = lue_daily_all(sub2);
% sify_daily_all2  = sify_daily_all(sub2);
% 
% figure
% scatter(lue_daily_all(sub1),sify_daily_all(sub1),24,apar_daily_all(sub1),'filled')
% plot(lue_daily_all(sub),sify_daily_all(sub),'r.','MarkerSize',30);
% lsline
% xlim([0,0.15])
% 
% xl          = xlim;
% yl          = ylim;
% textpos     = [0.05*(xl(2)-xl(1)) + xl(1), yl(2) - 0.05*(yl(2)-yl(1))];
% 
% fitresult   = fitlm(lue_daily_all(sub),sify_daily_all(sub));
% hh          = text(textpos(1),textpos(2),['r^{2}=' num2str(fitresult.Rsquared.Ordinary,'%4.2f')]);
% hh.FontSize = 16;
% 
% xlabel('LUE','FontSize',20);
% ylabel('SIFy','FontSize',20);
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/lue_sify_daily_scatter_all.png')
% 
% close(gcf)
% 
% 
% figure
% plot(lue_hourly_all,sify_hourly_all,'ro','MarkerSize',12);
% 
% lsline
% 
% xl          = xlim;
% yl          = ylim;
% textpos     = [0.05*(xl(2)-xl(1)) + xl(1), yl(2) - 0.05*(yl(2)-yl(1))];
% 
% fitresult   = fitlm(lue_hourly_all,sify_hourly_all);
% hh          = text(textpos(1),textpos(2),['r^{2}=' num2str(fitresult.Rsquared.Ordinary,'%4.2f')]);
% hh.FontSize = 16;
% 
% xlabel('LUE','FontSize',20);
% ylabel('SIFy','FontSize',20);
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/lue_sify_hourly_scatter_all.png')
% 
% close(gcf)

%% Fig. GPP EVI SIF
sif_daily_2014(sif_daily_2014>=1.0) = NaN;
[xData, yData] = prepareCurveData( doy_tmp_2014, sif_daily_2014 );

% Set up fittype and options.
ft = fittype( 'a+(b-k*x)*(1/(1+exp(c-d*x))-1/(1+exp(e-f*x)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0 0];
opts.StartPoint = [0.05 0.7 17.25 0.1128 33.1 0.117 0.005];
opts.Upper = [0.2 1 Inf Inf Inf Inf Inf];

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

[xData, yData] = prepareCurveData( doy_tmp_2014, gpp_daily_2014 );

% Set up fittype and options.
ft = fittype( 'a+(b-k*x)*(1/(1+exp(c-d*x))-1/(1+exp(e-f*x)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 15 150 1 270 1 0.005];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

[xData, yData] = prepareCurveData( doy_tmp_2014, evi_daily_2014 );

% Set up fittype and options.
ft = fittype( 'a+(b-k*x)*(1/(1+exp(c-d*x))-1/(1+exp(e-f*x)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0 0];
opts.StartPoint = [0.3 0.5 150 1 280 1 0.05];

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

%%% For Job Talk Figure %%%%%

font_size = 24;

% Step 1
hold on

figure('units','normalized','position',[0 0 1 1]);
h1     = gca;
set(h1,'Position',h1.Position - [0 0 0.1 0],...
       'YColor','k',...
       'XAxisLocation','bottom',... 
       'YAxisLocation','left',...
       'FontSize',font_size,...
       'FontName','Whitney',...
       'ylim',[0.2,0.7],...
       'xlim',[120,310],...
       'NextPlot','add');
plot(doy_tmp_2014,evi_daily_2014, 'k-');
ylabel('EVI','FontName','Whitney');
l1 = xlabel('Day of Year (2014)','FontSize',font_size,'FontName','Whitney');
h1_pos = h1.Position;
% a = get(gcf,'OuterPosition');
% set(gcf,'OuterPosition',[a(1),a(2),a(3)+0.5,a(4)+0.5*a(4)/a(3)]);
% set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,8]) % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/EVI_GPP_SIF_fig_step1.png');
%close(gcf);

% Step 2

xx3 = 120:1:310;
fit3 = fitresult{3};
yy3 = feval(fit3,xx3);
line(xx3,yy3,'Parent',h1,'LineWidth',2,'Color','k');
plot(141.4,0.2,'kv','MarkerSize',20,'Parent',h1,'MarkerFaceColor','k')
plot(291.2,0.2,'kv','MarkerSize',20,'Parent',h1,'MarkerFaceColor','k')

a = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[a(1),a(2),a(3)+0.5,a(4)+0.5*a(4)/a(3)]);
set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,8]) % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/EVI_GPP_SIF_fig_step2.png');
%close(gcf);


%'Position',get(h1,'Position')+[0 0 h1_pos(1,3) 0],...
h3     = axes( 'Position',h1_pos,...
                'XAxisLocation','top',...
               'XTick',[],...
               'XTickLabel',[],...
               'YAxisLocation','right',...                       
               'Color','none',...
               'XColor',get(gcf,'Color'),'YColor','b',...
               'FontSize',font_size,...
               'FontName','Whitney-Book',...
               'ylim',[0,25],...
               'xlim',[120,310]    );
ylabel('GPP (umol/m^{2}/second)','FontName','Whitney','Color','b');
[hb3,p3] = boundedline(doy_tmp_2014,gpp_daily_2014, zeros(size(gpp_daily_2014)), 'b-', h3,'alpha','transparency', 0);
%plot(doy_tmp_2014,gpp_daily_2014, 'k-');
plot(152.9,0,'bv','MarkerSize',20,'Parent',h3,'MarkerFaceColor','b')
plot(283.0,0,'bv','MarkerSize',20,'Parent',h3,'MarkerFaceColor','b')

xx2 = 120:1:300;
fit2 = fitresult{2};
yy2 = feval(fit2,xx2);
line(xx2,yy2,'Parent',h3,'LineWidth',2,'Color','b');

a = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[a(1),a(2),a(3)+0.5,a(4)+0.5*a(4)/a(3)]);
set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,8]) % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/EVI_GPP_SIF_fig_step3.png');
%close(gcf);
% get(h1,'Position')+[0 0 (20/190)*h1_pos(1,3) 0]
h2     = axes( 'Position',get(h1,'Position')+[0 0 (20/190)*h1_pos(1,3) 0],...
               'XAxisLocation','top',...
               'XTick',[],...
               'XTickLabel',[],...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','r',...
               'FontSize',font_size,...
               'FontName','Whitney-Book',...
               'ylim',[0,1],...
               'xlim',[120,330]    );
ylabel(h2, 'SIF(mw/m^{2}/nm/sr)','FontName','Whitney');
[hb2,p2] = boundedline(doy_tmp_2014,sif_daily_2014, zeros(size(sif_daily_2014)), 'r-', h2,'alpha','transparency', 0);
set(hb2,'MarkerSize',16);

plot(149.9,0,'ro','MarkerSize',20,'Parent',h2,'MarkerFaceColor','r')
plot(282.9,0,'ro','MarkerSize',20,'Parent',h2,'MarkerFaceColor','r')
xx1 = 120:1:300;
fit1 = fitresult{1};
yy1 = feval(fit1,xx1);
line(xx1,yy1,'Parent',h2,'LineWidth',2,'Color','r');

a = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[a(1),a(2),a(3)+0.5,a(4)+0.5*a(4)/a(3)]);
set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,8]) % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/EVI_GPP_SIF_fig_step4.png');
close(gcf);



hold off
        



%%% For Job Talk Figure %%%%%



% %%%%%--------------------%%%%%%%%%%%%%
% 
% % ax = subplot(2,1,1);
% % 
% % P = get(ax,'pos');    % Get the position.
% % delete(ax)            % Delete the subplot axes
% 
% FontName = 'Whitney'; %Gill Sans MT
% FontSize = 16;
% TickSize = 12;
% hold on
% %set(gca,'Color',[0 0 0])
% 
% [ax,h] = plotyyy(doy_tmp_2014,sif_daily_2014,doy_tmp_2014,gpp_daily_2014,doy_tmp_2014,evi_daily_2014);
% set(h(1),'Color','r')%,'marker','.','MarkerSize',20)
% %set(h,'linestyle','none','marker','.','MarkerSize',20)
% 
% set(h(2),'Color','k')%,'marker','.','MarkerSize',20);
% set(h(3),'Color','b')%,'marker','.','MarkerSize',20)
% 
% xlabel(ax(1),'Day of Year','FontName',FontName,'FontSize',FontSize);
% set(ax(1),'FontName',FontName,'FontSize',FontSize,'YColor','r');
% set(ax(2),'FontName',FontName,'FontSize',FontSize,'YColor','r');
% set(ax(3),'FontName',FontName,'FontSize',FontSize,'YColor','r');
% 
% %axes(ax(1,1))
% xx1 = 120:1:300;
% fit1 = fitresult{1};
% yy1 = feval(fit1,xx1);
% line(xx1,yy1,'Parent',ax(1),'LineWidth',2,'Color','r');
% ylabel(ax(1),'SIF(mw/m^{2}/nm/sr)','FontName',FontName,'FontSize',FontSize)
% 
% 
% 
% xx2 = 120:1:300;
% fit2 = fitresult{2};
% yy2 = feval(fit2,xx2);
% line(xx2,yy2,'Parent',ax(2),'LineWidth',2,'Color','k');
% ylabel(ax(2),'GPP(umol/m^{2}/sec)','FontName',FontName,'FontSize',FontSize)
% set(ax(2),'FontName',FontName,'FontSize',FontSize,'YColor','k');
% 
% 
% 
% xx3 = 120:1:300;
% fit3 = fitresult{3};
% yy3 = feval(fit3,xx3);
% line(xx3,yy3,'Parent',ax(3),'LineWidth',2,'Color','b');
% 
% ylabel(ax(3),'EVI','FontName',FontName,'FontSize',FontSize)
% set(ax(3),'FontName',FontName,'FontSize',FontSize,'YColor','b','YLim',[0.2,0.7]);
% 
% set(ax(1),'pos',P+[0 0 -0.7/5.5 0])
% set(ax(2),'pos',P+[0 0 -0.7/5.5 0])
% set(ax(3),'pos',P+[0 0 0 0])
% 
% plot(149.9,0.02,'ro','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','r')
% plot(152.9,0.02,'kv','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','k')
% plot(141.4,0.02,'bv','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','b')
% plot(282.9,0.02,'ro','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','r')
% plot(283.0,0.02,'kv','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','k')
% plot(291.2,0.02,'bv','MarkerSize',12,'Parent',ax(1),'MarkerFaceColor','b')
% 
% 
% hold off
% 
% % 
% % print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/2014_GPP_SIF_EVI.png')
% % 
% % close(gcf)
% 
% 
% %subplot(2,1,2)
% 
% 
% % FontName = 'Gill Sans MT';
% % FontSize = 16;
% % TickSize = 12;
% % [ax2,h2] = plotyy(doy_tmp_2014,rin_daily_2014,doy_tmp_2014,airT_daily_2014);
% % set(ax2,'FontName',FontName,'FontSize',FontSize);
% % 
% % xlabel(ax2(1),'Day of Year','FontName',FontName,'FontSize',FontSize,'Color','w');
% % ylabel(ax2(1),'shortwave radiation(w/m^{2})','FontName',FontName,'FontSize',FontSize)
% % ylabel(ax2(2),'air temperature(C)','FontName',FontName,'FontSize',FontSize)
% % 
% % set(ax2(1),'pos',P+[0 -0.45 -0.7/5.5 0])
% % set(ax2(2),'pos',P+[0 -0.45 -0.7/5.5 0])
% % a = get(gcf,'OuterPosition');
% % set(gcf,'OuterPosition',[a(1),a(2),a(3)+200,a(4)]);        %a(4)+0.2*a(4)/a(3)
% % set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,100,8])
% %set(gca,'Color',[0 0 0])
% export_fig('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/Phenology_EVI_SIF_GPP.png','-r300','-dpng','-transparent','-nocrop')
% 
% %print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/2014_airT_Rin_EVI_SIF_GPP.png')
% 
% close(gcf)

















