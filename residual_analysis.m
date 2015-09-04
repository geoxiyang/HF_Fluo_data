%% GPP = f(APAR) + residual1; SIF = f(APAR) + residual2. resid1 ~ resid2


datapath = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';
load([datapath,'SIF760daily.mat'],'halfhourly_result','raw_final_result','SIF_mean');
load([datapath,'HF_2013_GPP.mat']); 
load([datapath,'hf_barn_2013_env.mat'],'apar','apar_daily')

%% Plot

% Fontsize = 16;
% opath    = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/';
% 
% % GPP
% figure
% plot(apar(gpp_raw>0),gpp_raw(gpp_raw>0),'ko')
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('APAR(umol/m^{2}/second)','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('GPP(umol/m^{2}/second)','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'HF_GPP_APAR_30min.png']);
% close(gcf);
% 
% 
% figure
% plot(apar_daily,gpp_day,'ko')
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('APAR(umol/m^{2}/second)','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('GPP(g C/m^{2}/day)','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'HF_GPP_APAR_daily.png']);
% close(gcf);
% 
% % SIF
% figure
% plot(apar(halfhourly_result(:,2)>0),halfhourly_result(halfhourly_result(:,2)>0,2),'ro')
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('APAR(umol/m^{2}/second)','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('SIF(mw/m^{2}/sr/nm)','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'HF_SIF_APAR_30min.png']);
% close(gcf);
% 
% 
% figure
% plot(apar_daily(SIF_mean(:,2)>0),SIF_mean(SIF_mean(:,2)>0,2),'ro')
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('APAR(umol/m^{2}/second)','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('SIF(mw/m^{2}/sr/nm)','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'HF_SIF_APAR_daily.png']);
% close(gcf);

%% Residual analysis

Fontsize = 16;
opath    = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/';

sub     = gpp_raw>0 & halfhourly_result(:,2)>0 & ~isnan(apar);

f1      = polyfit(apar(sub),gpp_raw(sub),1);
yfit1   = polyval(f1,apar(sub));
res1    = gpp_raw(sub) - yfit1;

f2      = polyfit(apar(sub),halfhourly_result(sub,2),1);
yfit2   = polyval(f2,apar(sub));
res2    = halfhourly_result(sub,2) - yfit2;

% plot(res1,res2,'ko')
% 
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('GPP residual','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('SIF residual','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'GPP_APAR_30min_residual_linear.png']);
% close(gcf);

sub2    = gpp_day>0 & SIF_mean(:,2)>0 & ~isnan(apar_daily);

f3      = polyfit(apar_daily(sub2),gpp_day(sub2),1);
yfit3   = polyval(f3,apar_daily(sub2));
res3    = gpp_day(sub2) - yfit3;

f4      = polyfit(apar_daily(sub2),SIF_mean(sub2,2),1);
yfit4   = polyval(f4,apar_daily(sub2));
res4    = SIF_mean(sub2,2) - yfit4;

% plot(res3,res4,'ko')
% set(gca,'FontSize',Fontsize,'FontName','Whitney')
% xlabel('GPP residual','FontName','Whitney','Color','k','FontSize',Fontsize)
% ylabel('SIF residual','FontName','Whitney','Color','k','FontSize',Fontsize)
% print(gcf, '-dpng','-r300', [opath 'GPP_APAR_daily_residual_linear.png']);
% close(gcf);









