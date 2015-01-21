%% Calculate daily LUE and SIFyield

% Note: the method I used in previous versions of the manuscript is that
% I calculated LUE every 30min and used only the mean of midday LUE as the
% daily LUE.
% SIFyield is daily SIF/daily APAR

datapath = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';

load([datapath,'hf_barn_2013_env.mat'],'apar_daily','apar','doy','par','cloud_ratio','cloud_ratio_daily')
load([datapath,'HF_2013_GPP.mat'],'gpp_day','gpp_raw')
load([datapath,'SIF760daily.mat'],'SIF_mean','SIF_1330','SIF_0930','SIF_1400','halfhourly_result')

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

% 30min LUE and SIFyield
day_sub_30min       = find(apar >=100);
sif_yield_30min     = nan(numel(apar),1);
LUE_30min           = nan(numel(apar),1);
sif_yield_30min(day_sub_30min) = halfhourly_result(day_sub_30min,2)./apar(day_sub_30min);
LUE_30min(day_sub_30min) = gpp_raw(day_sub_30min)./apar(day_sub_30min);

% SIFyield and LUE should be larger than 0
sif_yield_30min(sif_yield_30min<=0) = NaN;
LUE_30min(LUE_30min<=0) = NaN;

% daily LUE and SIFyield
LUE_day  = nan(130,3);
SIFy_day = nan(130,3);
for jj = 1:130
    LUE_day(jj,1) = nanmean(LUE_30min(doy>= 169+jj & doy < 169 + jj + 1));
    SIFy_day(jj,1)= nanmean(sif_yield_30min(doy>= 169+jj & doy < 169 + jj + 1));
    if cloud_ratio_daily(jj,1) < 0.5
        % less than 50% diffuse light        
        LUE_day(jj,2) = nanmean(LUE_30min(doy>= 169+jj & doy < 169 + jj + 1));
        SIFy_day(jj,2)= nanmean(sif_yield_30min(doy>= 169+jj & doy < 169 + jj + 1));
        LUE_day(jj,3) = NaN;
        SIFy_day(jj,3)= NaN;
    else
        % more than 50% diffuse light
        LUE_day(jj,2) = NaN;
        SIFy_day(jj,2)= NaN;
        LUE_day(jj,3) = nanmean(LUE_30min(doy>= 169+jj & doy < 169 + jj + 1));
        SIFy_day(jj,3)= nanmean(sif_yield_30min(doy>= 169+jj & doy < 169 + jj + 1));
    end
end

% diffuse proportion vs. LUE_day & SIFy_day
%figure('units','normalized','position',[0 0 1 1])

% [AX,H1,H2] = plotyy(cloud_ratio_daily,LUE_day(:,1),cloud_ratio_daily,SIFy_day(:,1));
% ylabel(AX(1),'LUE(umol CO_{2}/umol photon)','FontName','Whitney','FontSize',20,'Color','r');
% ylabel(AX(2),'SIF/APAR','FontName','Whitney','FontSize',20,'Color','m');
% set(AX(1),'YColor','r','FontSize',16,'FontName','Whitney');
% set(AX(2),'YColor','m','FontSize',16,'FontName','Whitney');
% set(H1,'linestyle','none','marker','o','color','r','MarkerSize',12);
% set(H2,'linestyle','none','marker','^','color','m','MarkerSize',12);
% xlabel(AX(1),'Diffuse PAR Fraction','FontName','Whitney','FontSize',20);
% set(gcf,'paperPositionMode','auto')
% print(gcf, '-dpng','-r300', '/Users/xiyang/Dropbox/Mypaper/6.2014-Fluorescence/resub2/figures/DiffuseFrac_SIFy_LUE.png')


% Daily map ========%
subs = LUE_day > 0 & SIFy_day > 0;
LUE_day(~subs)   = NaN;
SIFy_day(~subs)  = NaN;

subs2 = LUE_day(:,1) > 0 & SIFy_day(:,1) >0;
subs3 = LUE_day(:,2) > 0 & SIFy_day(:,2) >0;
subs4 = LUE_day(:,3) > 0 & SIFy_day(:,3) >0;
% 
% APAR vs LUE and SIFy
% [AX,H1,H2] = plotyy(apar_daily(subs2),LUE_day(subs2,1),apar_daily(subs2),SIFy_day(subs2,1));
% ylabel(AX(1),'LUE(umol CO_{2}/umol photon)','FontName','Whitney','FontSize',20);
% ylabel(AX(2),'SIF/APAR','FontName','Whitney','FontSize',20);
% set(H1,'linestyle','none','marker','o','color','r');
% set(H2,'linestyle','none','marker','o','color','m');
% xlabel(AX(1),'APAR(umol/m^{2}/second)','FontName','Whitney','FontSize',20);
% set(gcf,'paperPositionMode','auto')
% 
% hold on 
% [AX2,H12,H22] = plotyy(apar_daily(subs4),LUE_day(subs4,1),apar_daily(subs4),SIFy_day(subs4,1));
% ylabel(AX2(1),'LUE(umol CO_{2}/umol photon)','FontName','Whitney','FontSize',20);
% ylabel(AX2(2),'SIF/APAR','FontName','Whitney','FontSize',20);
% set(H12,'linestyle','none','marker','.','color','r');
% set(H22,'linestyle','none','marker','.','color','m');
% xlabel(AX2(1),'APAR(umol/m^{2}/second)','FontName','Whitney','FontSize',20);
% 
% 
% x =  [ones(numel(apar_daily(sub2,1)),1),LUE_day(subs2,1)];
% [b1,bint1,r1,rint1,stats1] = regress(LUE_day(subs2,1),x);
% [b2,bint2,r2,rint2,stats2] = regress(SIFy_day(subs2,1),x);
% 
% xint1 = (max(LUE_day(subs2,1)) - min(LUE_day(subs2,1)))/100.0;
% xfit1 = min(LUE_day(subs2,1)):xint:max(LUE_day(subs2,1));
% yfit1 = b(1) + b(2) * xfit;
% 
% 
% 

% indice = 1:125;
% plot(LUE_day(indice,1),SIFy_day(indice,1),'ko','MarkerSize',16)
% corr(LUE_day(indice,1),SIFy_day(indice,1))^2
% x = [ones(numel(LUE_day(indice,1)),1),LUE_day(indice,1)];
% [b,bint,r,rint,stats] = regress(SIFy_day(indice,1),x);
% 
% xlabel('LUE(umol CO_{2}/umol photon)','FontName','Whitney','FontSize',20)
% ylabel('SIF/APAR','FontName','Whitney','FontSize',20)
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') 
% 
% xlim_fig = xlim;
% ylim_fig = ylim;
% 
% wdth = xlim_fig(2)-xlim_fig(1);
% ht = ylim_fig(2)-ylim_fig(1);
% pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
% text(pos(1),pos(2),'R^{2} = 0.3868     p=0.0000','FontSize',20,'FontName','Whitney-Book');
% 
% hold on
% xint = (max(LUE_day(indice,1)) - min(LUE_day(indice,1)))/100.0;
% xfit = min(LUE_day(indice,1)):xint:max(LUE_day(indice,1));
% yfit = b(1) + b(2) * xfit;
% 
% plot(xfit,yfit,'b-');
% 
% hold off

% 
% plot(LUE_day(subs2),SIFy_day(subs2),'ko','MarkerSize',16)
% corr(LUE_day(subs2),SIFy_day(subs2))^2
% x = [ones(numel(LUE_day(subs2,1)),1),LUE_day(subs2,1)];
% [b,bint,r,rint,stats] = regress(SIFy_day(subs2,1),x);
% 
% xlabel('LUE(umol CO_{2}/umol photon)','FontName','Whitney','FontSize',20)
% ylabel('SIF/APAR','FontName','Whitney','FontSize',20)
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') 
% 
% xlim_fig = xlim;
% ylim_fig = ylim;
% 
% wdth = xlim_fig(2)-xlim_fig(1);
% ht = ylim_fig(2)-ylim_fig(1);
% pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
% text(pos(1),pos(2),'R^{2} = 0.3868     p=0.0000','FontSize',20,'FontName','Whitney-Book');
% 
% hold on
% xint = (max(LUE_day(subs2,1)) - min(LUE_day(subs2,1)))/100.0;
% xfit = min(LUE_day(subs2,1)):xint:max(LUE_day(subs2,1));
% yfit = b(1) + b(2) * xfit;
% 
% plot(xfit,yfit,'b-');
% 
% %======================%
% 
% % sub_3 = cloud_ratio < 0.5  & LUE_30min > 0 & sif_yield_30min > 0 & apar > 1000;  %&
% % scatter(LUE_30min(sub_3),sif_yield_30min(sub_3),[],apar(sub_3))
% 
% 
% 
% 
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIFy_LUE.png')

% LUE      = gpp_day./apar_daily;
% SIFyield = SIF_mean(:,2)./apar_daily;
% 
% SIFyield_sun = SIF_mean(:,3)./apar_daily;
% SIFyield_cloud = SIF_mean(:,4)./apar_daily;
% 
% save('SIFyield_2013_SVD.mat');

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


