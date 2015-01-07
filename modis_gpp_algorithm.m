%% Calculate daily GPP using MODIS algorithm; Plot comparison between MODIS GPP (algorithm & satellite) vs. EC GPP
% MOD17a1
% GPP (g C m-2) = LUE_max * 1e3 * Rshort * 86400 * 1e-6 * 0.45 * fAPAR * Tmin_SCALAR * VPD_SCALAR

clear all

load('HF_2013_GPP.mat','gpp_day');

%GPP Data from MODIS
%Integrate 8-day GPP EC data to match MODIS GPP

% sat_gpp = importdata('/Volumes/XiYangResearch/Projects/9.Fluorescence/8.SatelliteData/modis/gpp_hf_2013.txt');
% %sat_qc  = importdata('/Volumes/XiYangResearch/Projects/9.Fluorescence/8.SatelliteData/modis/gpp_qc_hf_2013.txt');
% 
% sat_gpp(:,2) = sat_gpp(:,2) * 0.0001 * 1000;
% sat_gpp(:,1) = (sat_gpp(:,1)-1)*8+1;
% for ii = 1:numel(sat_gpp(:,1))
%    date_range = [(ii-1)*8+1,ii*8];
%    if date_range(1) >= 170 && date_range(2) <= 300
%       ec_gpp_8day(ii) = sum(gpp_day((date_range(1)-170+1):(date_range(2)-170+1)));
%    else
%       ec_gpp_8day(ii) = NaN;
%    end
% end
% 
% % %1. Plot Time-series of MODIS & EC GPP
% plot(sat_gpp(:,1),sat_gpp(:,2),'ro',sat_gpp(:,1),ec_gpp_8day,'bo','MarkerSize',16);
% xlim([170,300]);
% ylim([0,100]);
% xlabel('Day of Year','FontName','Whitney','FontSize',24);
% ylabel('8-day MODIS GPP & EC GPP (g C/m^{2})','FontName','Whitney','FontSize',24);
% set(gca,'FontName','Whitney','FontSize',20);
% legend('MODIS GPP','EC GPP');
% 
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/SAT_EC_GPP_TS.png');
% close(gcf);

% %2. Plot scatter plot of MODIS & EC GPP
% plot(ec_gpp_8day,sat_gpp(:,2),'ko','MarkerSize',16)
% hold on
% xint = (nanmax(ec_gpp_8day) - nanmin(ec_gpp_8day))/100.0;
% xfit = nanmin(ec_gpp_8day):xint:nanmax(ec_gpp_8day);
% x    = [ones(numel(ec_gpp_8day(~isnan(ec_gpp_8day))),1),ec_gpp_8day(~isnan(ec_gpp_8day))'];
% [b,bint,r,rint,stats] = regress(sat_gpp(~isnan(ec_gpp_8day),2),x);
% yfit = b(1) + b(2) * xfit;
% plot(xfit,yfit,'b-');
% 
% xlabel('8-day EC GPP(g C/m^{2})','FontName','Whitney','FontSize',24);
% ylabel('8-day MODIS GPP(g C/m^{2})','FontName','Whitney','FontSize',24);
% set(gca,'FontName','Whitney','FontSize',20);
% title('');
% 
% pos = [2,75];
% text(pos(1),pos(2),'r^{2}=0.502','FontSize',24,'FontName','Whitney');
% hold off
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/SAT_EC_GPP_SC.png');
% close(gcf);


% GPP Algorithm & local meterological data
load('hf_barn_2013_env.mat','vpd_daily','tmin_daily','Rshort_daily','fapar_daily','apar_daily')
LUE_max = 0.001044;

Tmin_SCALAR = zeros(130,1);
VPD_SCALAR  = zeros(130,1);

for ii = 1:130
    if tmin_daily(ii) < -8.00
        Tmin_SCALAR(ii) = 0;
    elseif tmin_daily(ii) > 7.94
        Tmin_SCALAR(ii) = 1;
    else
        Tmin_SCALAR(ii) = (tmin_daily(ii) - (-8.00))/(7.94 - (-8.00));
    end
    
    if vpd_daily(ii) > 3100
        VPD_SCALAR(ii) = 0;
    elseif vpd_daily(ii) < 650
        VPD_SCALAR(ii) = 1;
    else
        VPD_SCALAR(ii) = (3100 - vpd_daily(ii))/(3100 - 650);
    end   
    
end

GPP_MODIS_ALG = LUE_max * 1e3 * 1e-6 .* Rshort_daily * 86400 * 0.45 .* fapar_daily .* Tmin_SCALAR .* VPD_SCALAR;

%GPP_MODIS_ALG2 = LUE_max * 1e3 * 1e-6 .* apar_daily * 86400 * 0.219 .* fapar_daily .* Tmin_SCALAR .* VPD_SCALAR;

GPP_MODIS_ALG(GPP_MODIS_ALG >40)    = NaN;
%GPP_MODIS_ALG2(GPP_MODIS_ALG2 >40)  = NaN;

% %3. Plot scatter plot of MODIS & EC ALGORITHM GPP
plot(gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1),GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1),'ko','MarkerSize',16)
hold on

xint = (nanmax(gpp_day(:,1)) - nanmin(gpp_day(:,1)))/1000.0;
xfit = nanmin(gpp_day(:,1)):xint:nanmax(gpp_day(:,1));
x    = [ones(numel(gpp_day(~isnan(GPP_MODIS_ALG(:,1)))),1),gpp_day(~isnan(GPP_MODIS_ALG(:,1)))];
[b,bint,r,rint,stats] = regress(GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1),x);
yfit = b(1) + b(2) * xfit;
plot(xfit,yfit,'b-');

xlabel('Daily EC GPP(g C/m^{2})','FontName','Whitney','FontSize',24);
ylabel('Daily MODIS GPP(g C/m^{2})','FontName','Whitney','FontSize',24);
set(gca,'FontName','Whitney','FontSize',20);
title('');

pos = [0.5,14];
text(pos(1),pos(2),'r^{2}=0.612','FontSize',24,'FontName','Whitney');
hold off

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/Fig.S4.png');
close(gcf);

%4. Plot time-series of MODIS algorithm GPP vs. EC GPP
% plot(170:1:299,GPP_MODIS_ALG(:,1),'ro',170:1:299,gpp_day(:,1),'bo','MarkerSize',16);
% xlim([170,300]);
% ylim([0,20]);
% xlabel('Day of Year','FontName','Whitney','FontSize',24);
% ylabel('Daily MODIS GPP & EC GPP (g C/m^{2})','FontName','Whitney','FontSize',24);
% set(gca,'FontName','Whitney','FontSize',20);
% legend('MODIS GPP','EC GPP');
% 
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/ALG_EC_GPP_TS.png');
% close(gcf)



% plot(gpp_day(~isnan(GPP_MODIS_ALG2(:,1)),1),GPP_MODIS_ALG2(~isnan(GPP_MODIS_ALG2(:,1)),1),'ro')
% hold off

% 
% aa = gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1);
% bb = GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1);
% 
% corr(gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1),GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1))^2

% plot(170:1:299,gpp_day(:,1),'ro',170:1:299,GPP_MODIS_ALG(:,1),'bo')