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
             'vpd_hourly');
         
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
             'vpd_hourly');
         
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
             'vpd_hourly');
%          
% need to make daily SIF, GPP, LE, APAR, PRI. Here we use data of all kinds of
% sky conditions

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


for ii = 1:length(doy_tmp_2012)
    
    gpp_daily_2012(ii,1)    = nanmean(data2012.GPP(data2012.GPP >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<(doy_tmp_2012(ii)+1)));
    le_daily_2012(ii,1)     = nanmean(data2012.LE(data2012.LE >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    sif_daily_2012(ii,1)    = nanmean(data2012.sif_hourly(data2012.sif_hourly >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    apar_daily_2012(ii,1)   = nanmean(data2012.apar_hourly(data2012.apar_hourly>0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    vpd_daily_2012(ii,1)    = nanmean(data2012.vpd_hourly(data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    
end

for ii = 1:length(doy_tmp_2013)
    
    gpp_daily_2013(ii,1)    = nanmean(data2013.GPP(data2013.GPP >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    le_daily_2013(ii,1)     = nanmean(data2013.LE(data2013.LE >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    sif_daily_2013(ii,1)    = nanmean(data2013.sif_hourly(data2013.sif_hourly >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    apar_daily_2013(ii,1)   = nanmean(data2013.apar_hourly(data2013.apar_hourly>0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    vpd_daily_2013(ii,1)    = nanmean(data2013.vpd_hourly(data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    
end

for ii = 1:length(doy_tmp_2014)
    
    gpp_daily_2014(ii,1)    = nanmean(data2014.GPP(data2014.GPP >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    le_daily_2014(ii,1)     = nanmean(data2014.LE(data2014.LE >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    sif_daily_2014(ii,1)    = nanmean(data2014.sif_hourly(data2014.sif_hourly >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    apar_daily_2014(ii,1)   = nanmean(data2014.apar_hourly(data2014.apar_hourly>0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    vpd_daily_2014(ii,1)    = nanmean(data2014.vpd_hourly(data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    
end

gpp_daily_all   = [gpp_daily_2012;gpp_daily_2013;gpp_daily_2014];
le_daily_all    = [le_daily_2012;le_daily_2013;le_daily_2014];
apar_daily_all  = [apar_daily_2012;apar_daily_2013;apar_daily_2014];
sif_daily_all   = [sif_daily_2012;sif_daily_2013;sif_daily_2014];
doy_daily_all   = [doy_tmp_2012;doy_tmp_2013+366;doy_tmp_2014+366+365];
vpd_daily_all   = [vpd_daily_2012;vpd_daily_2013;vpd_daily_2014];


gpp_hourly_all   = [data2012.GPP;data2013.GPP;data2014.GPP];
le_hourly_all    = [data2012.LE;data2013.LE;data2014.LE];
apar_hourly_all  = [data2012.apar_hourly;data2013.apar_hourly;data2014.apar_hourly];
sif_hourly_all   = [data2012.sif_hourly;data2013.sif_hourly;data2014.sif_hourly];
doy_hourly_all   = [data2012.doy_hourly;data2013.doy_hourly+366;data2014.doy_hourly+366+365];
sp_hourly_all    = [data2012.sunportion_hourly;data2013.sunportion_hourly;data2014.sunportion_hourly];

th = 0;

gpp_hourly_all(gpp_hourly_all<=0 | sp_hourly_all<th)     = NaN;
le_hourly_all(le_hourly_all<=0 | sp_hourly_all<th)       = NaN;
apar_hourly_all(apar_hourly_all<=0 | sp_hourly_all<th)   = NaN;
sif_hourly_all(sif_hourly_all<=0 | sp_hourly_all<th)     = NaN;

%% New code making Fig.1

% fontsize = 12;
% 
% subplot(4,1,1)
% plot(doy_daily_all,gpp_daily_all,'ro')
% % hold on
% % line([365,365],[0,30],'Color','k','LineStyle','--')
% % line([731,731],[0,30],'Color','k','LineStyle','--')
% % text(290,25,num2str(2012))
% % text(655,25,num2str(2013))
% % text(1022,25,num2str(2012))
% 
% % hold off
% set(gca,'FontSize',fontsize)
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel('GPP(umol/m^{2}/sec)','FontSize',fontsize);
% 
% subplot(4,1,2)
% plot(doy_daily_all,le_daily_all,'bo')
% set(gca,'FontSize',fontsize,'YLim',[0, 200])
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel('LE(w/m^{2})','FontSize',fontsize);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% subplot(4,1,3)
% plot(doy_daily_all,sif_daily_all,'mo')
% set(gca,'FontSize',fontsize)
% %xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel('SIF(mw/m^{2}/nm/sr)','FontSize',fontsize);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% subplot(4,1,4)
% plot(doy_daily_all,apar_daily_all,'ko')
% set(gca,'FontSize',fontsize)
% xlabel('Day of Year since 2012-01-01','FontSize',fontsize);
% ylabel('APAR(umol/m^{2}/sec)','FontSize',fontsize);%LE(w/m^{2})%GPP(umol/m^{2}/sec)  
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/Fig.1_1.png')
% close(gcf)


%% Fig.2 

% %---------------Fig.2a--------------------
% %-----------------------------------------
% figure
% 
% 
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
% h = plot(fitresult);
% set(h,'Color','k','LineWidth',2)
% set(gca,'YLim',[0,max(gpp_daily_all)*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,max(gpp_daily_all)*1.1,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
% hh.FontSize = 20;   
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
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,max(gpp_daily_all)*1.0,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('GPP(umol/m^{2}/sec)','FontSize',20);
% hold off
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_SIF_GPP.png')
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
% plot(sif_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all)),apar_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all)),'m.','MarkerSize',20);
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
% testa = sif_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all));
% testb = apar_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all));
% % Fit model to data.
% [fitresult3, gof3] = fit(testa, testb, ft3, opts3);
% 
% h = plot(fitresult3);
% set(h,'Color','k','LineWidth',2,'LineStyle','--');
% set(gca,'YLim',[0,max(apar_daily_all(~isnan(sif_daily_all) & ~isnan(apar_daily_all)))*1.2],...
%     'FontSize',16);
% legend('off')
% 
% hh = text(0.05,850,['Linear fit r^{2}=' num2str(gof3.rsquare,'%4.2f')]);
% hh.FontSize = 20; 
% 
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('APAR(umol/m^{2}/sec)','FontSize',20);
% hold off
% 
% print('-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_SIF_APAR.png')
% close(gcf)

%% Fig.2d-f

% 


filename1        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_GPP.png';
filename2        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_LE.png';
filename3        = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/scatter_annual_hourly_SIF_APAR.png';

%---------------Fig.2d--------------------
%-----------------------------------------
figure

% temp_hours = [data2012.doy_hourly-fix(data2012.doy_hourly);data2013.doy_hourly-fix(data2013.doy_hourly);data2014.doy_hourly-fix(data2014.doy_hourly)];
% temp_subs  = temp_hours>=9.9/24.0 & temp_hours<=10.1/24.0;
% 
% sif_hourly_all = sif_hourly_all(temp_subs);
% gpp_hourly_all = gpp_hourly_all(temp_subs);
% 
plot(sif_hourly_all,gpp_hourly_all,'r.','MarkerSize',20);




% Label axes
%LE(w/m^{2})%GPP(umol/m^{2}/sec)

hold on
% Michaelis-Menten fit
ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Off';
opts.StartPoint = [max(gpp_hourly_all) 0.5];
opts.Lower  = [0 0];
opts.Upper  = [2*max(gpp_hourly_all) 5];

% Fit model to data.
[fitresult, gof] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), gpp_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), ft, opts );

h = plot(fitresult);
set(h,'Color','r','LineWidth',2)
set(gca,'YLim',[0,max(gpp_hourly_all)*1.2],...
    'FontSize',16);
legend('off')

hh = text(0.05,max(gpp_hourly_all)*1.1,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
hh.FontSize = 20;   

% Linear fit
ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Off';
opts2.StartPoint = [20 1];
opts2.Lower  = [0 0];
opts2.Upper  = [2*max(gpp_hourly_all) 5];

% Fit model to data.
[fitresult2, gof2] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), gpp_hourly_all(~isnan(sif_hourly_all) & ~isnan(gpp_hourly_all)), ft2, opts2);

h = plot(fitresult2);
set(h,'Color','r','LineWidth',2,'LineStyle','--');
set(gca,'YLim',[0,max(gpp_hourly_all)*1.2],...
    'FontSize',16);
legend('off')

hh = text(0.05,max(gpp_hourly_all)*1.0,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
hh.FontSize = 20; 

xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
ylabel('GPP(umol/m^{2}/sec)','FontSize',20);
hold off

print('-dpng','-r300',filename1)
close(gcf)
%-----------------------------------------
%---------------Fig.2b--------------------
%-----------------------------------------
figure
plot(sif_hourly_all,le_hourly_all,'b.','MarkerSize',20);
% Label axes
%LE(w/m^{2})%GPP(umol/m^{2}/sec)

hold on
% Michaelis-Menten fit
ft = fittype('a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Off';
opts.StartPoint = [max(le_hourly_all) 0.5];
opts.Lower  = [0 0];
opts.Upper  = [2*max(le_hourly_all) 5];

% Fit model to data.
[fitresult, gof] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), le_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), ft, opts );

h = plot(fitresult);
set(h,'Color','k','LineWidth',2)
set(gca,'YLim',[0,max(le_hourly_all)*1.05],...
    'FontSize',16);
legend('off')

hh = text(0.05,500,['MM fit r^{2}=' num2str(gof.rsquare,'%4.2f')]);
hh.FontSize = 20;   

% Linear fit
ft2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Off';
opts2.StartPoint = [20 1];
opts2.Lower  = [0 0];
opts2.Upper  = [2*max(le_hourly_all) 5];

% Fit model to data.
[fitresult2, gof2] = fit( sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), le_hourly_all(~isnan(sif_hourly_all) & ~isnan(le_hourly_all)), ft2, opts2);

h = plot(fitresult2);
set(h,'Color','k','LineWidth',2,'LineStyle','--');
set(gca,'YLim',[0,max(le_hourly_all)*1.2],...
    'FontSize',16);
legend('off')

hh = text(0.05,460,['Linear fit r^{2}=' num2str(gof2.rsquare,'%4.2f')]);
hh.FontSize = 20; 

xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
ylabel('LE(w/m^{2})','FontSize',20);

set(gca,'YLim',[0,max(le_hourly_all)*1.05],...
    'FontSize',16);
hold off

print('-dpng','-r300',filename2)
close(gcf)

%-----------------------------------------
%---------------Fig.2c--------------------
%-----------------------------------------
figure
plot(sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)),apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)),'m.','MarkerSize',20);

hold on

% Linear fit
ft3 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Robust = 'Off';
opts3.StartPoint = [500 120];
opts3.Lower  = [-Inf -Inf];
opts3.Upper  = [Inf Inf];

testa = sif_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all));
testb = apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all));
% Fit model to data.
[fitresult3, gof3] = fit(testa, testb, ft3, opts3);

h = plot(fitresult3);
set(h,'Color','k','LineWidth',2,'LineStyle','--');
set(gca,'YLim',[0,max(apar_hourly_all(~isnan(sif_hourly_all) & ~isnan(apar_hourly_all)))*1.2],...
    'FontSize',16);
legend('off')

hh = text(0.05,2200,['Linear fit r^{2}=' num2str(gof3.rsquare,'%4.2f')]);
hh.FontSize = 20; 

xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
ylabel('APAR(umol/m^{2}/sec)','FontSize',20);
hold off

print('-dpng','-r300',filename3)
close(gcf)

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

% lue_daily_all   = gpp_daily_all./apar_daily_all;
% sify_daily_all  = sif_daily_all./apar_daily_all;
% lue_hourly_all  = gpp_hourly_all./apar_hourly_all;
% sify_hourly_all = sif_hourly_all./apar_hourly_all;
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



