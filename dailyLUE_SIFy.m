%% Calculate daily LUE and SIFyield

% Note: the method I used in previous versions of the manuscript is that
% I calculated LUE every 30min and used only the mean of midday LUE as the
% daily LUE.
% SIFyield is daily SIF/daily APAR


load('hf_barn_2013_env.mat','apar_daily','apar','doy','par')
load('HF_2013_GPP.mat','gpp_day')
load('SIF760daily_2013_SVD.mat','SIF_mean','SIF_1330','SIF_0930','SIF_1400')

apar_1330 = zeros(130,1);
apar_1400 = zeros(130,1);
apar_0930 = zeros(130,1);
par_0930  = zeros(130,1);

for uni_i = 1:130
   
   apar_1330(uni_i,1) = apar(doy>=uni_i+169+13.4/24.0 & doy<uni_i+169+13.6/24.0);
   apar_0930(uni_i,1) = apar(doy>=uni_i+169+09.4/24.0 & doy<uni_i+169+09.6/24.0);
   apar_1400(uni_i,1) = apar(doy>=uni_i+169+13.9/24.0 & doy<uni_i+169+14.1/24.0);
   par_0930(uni_i,1)  = par(doy>=uni_i+169+09.4/24.0 & doy<uni_i+169+09.6/24.0);
   
end



sif_yield_1330(:,1) = SIF_1330(:,1);
sif_yield_1330(:,2) = SIF_1330(:,2)./apar_1330;
sif_yield_1330(:,3) = SIF_1330(:,3)./apar_1330;
sif_yield_1330(:,4) = SIF_1330(:,4)./apar_1330;

sif_yield_1400(:,1) = SIF_1400(:,1);
sif_yield_1400(:,2) = SIF_1400(:,2)./apar_1400;
sif_yield_1400(:,3) = SIF_1400(:,3)./apar_1400;
sif_yield_1400(:,4) = SIF_1400(:,4)./apar_1400;

sif_yield_0930(:,1) = SIF_0930(:,1);
sif_yield_0930(:,2) = SIF_0930(:,2)./apar_0930;
sif_yield_0930(:,3) = SIF_0930(:,3)./apar_0930;
sif_yield_0930(:,4) = SIF_0930(:,4)./apar_0930;


LUE      = gpp_day./apar_daily;
SIFyield = SIF_mean(:,2)./apar_daily;

SIFyield_sun = SIF_mean(:,3)./apar_daily;
SIFyield_cloud = SIF_mean(:,4)./apar_daily;

save('SIFyield_2013_SVD.mat');

%plot(sif_yield_0930(sif_yield_0930(:,2)>0,1),sif_yield_0930(sif_yield_0930(:,2)>0,2),'ro')

% plot(sif_yield_0930(sif_yield_0930(:,3)>0,1),sif_yield_0930(sif_yield_0930(:,3)>0,3),'ro','MarkerSize',12)
% hold on
% plot(sif_yield_0930(sif_yield_0930(:,4)>0,1),sif_yield_0930(sif_yield_0930(:,4)>0,4),'bo','MarkerSize',12)
% xlabel('Day of Year','FontName','Whitney','FontSize',20)
% ylabel('SIF/APAR','FontName','Whitney','FontSize',20)
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') 
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIFy_Season.png')


% plot(sif_yield_0930(sif_yield_0930(:,2)>0,2),VPD_0930(sif_yield_0930(:,2)>0,2),'ko');
% plot(sif_yield_1400(sif_yield_1400(:,2)>0,2),VPD_1400(sif_yield_1400(:,2)>0,2),'ko');


% figure
% plot(SIF_mean(SIF_mean(:,2)>0,2),gpp_day(SIF_mean(:,2)>0),'ko');
% 
% 
% figure
% plot(SIFyield(SIFyield>0),LUE(SIFyield>0),'bo')
% 
% hold on
% plot(LUE(SIFyield_sun>0),SIFyield_sun(SIFyield_sun>0),'ro','MarkerSize',12)
% plot(LUE(SIFyield_cloud>0),SIFyield_cloud(SIFyield_cloud>0),'bo','MarkerSize',12)
% 
% 
% corr(LUE(SIFyield_sun>0),SIFyield_sun(SIFyield_sun>0))^2
% corr(LUE(SIFyield_cloud>0),SIFyield_cloud(SIFyield_cloud>0))^2
% 
% hold off
% 
% xlabel('LUE','FontName','Whitney','FontSize',20);
% ylabel('SIF_{yield}','FontName','Whitney','FontSize',20);
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIFy_LUE.png');


