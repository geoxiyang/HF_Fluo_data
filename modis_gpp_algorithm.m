%% Calculate daily GPP using MODIS algorithm
% MOD17a1
% GPP (g C m-2) = LUE_max * 1e3 * Rshort * 86400 * 1e-6 * 0.45 * fAPAR * Tmin_SCALAR * VPD_SCALAR

clear all

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


load('HF_2013_GPP.mat','gpp_day');
plot(gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1),GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1),'ko')
hold on
h   = plotregression(gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1),GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1),'regression');
h1  = get(h,'Children');
h2  = get(h1(2),'Children');
delete(h2(3));     % delete y=T fit
delete(h1(1));     % delete legend box

xlabel('GPP EC','FontName','Whitney','FontSize',30);
ylabel(h1(2),'GPP MODIS','FontName','Whitney','FontSize',30);
set(h1(2),'FontName','Whitney','FontSize',24);
title('');

pos = [1,13];
text(pos(1),pos(2),'r^{2}=0.612','FontSize',30,'FontName','Whitney');
hold off

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Users/xiyang/Dropbox/Mypaper/6.2014-Fluorescence/JPG/Final/Fig.S4.png');
close(gcf);


% plot(gpp_day(~isnan(GPP_MODIS_ALG2(:,1)),1),GPP_MODIS_ALG2(~isnan(GPP_MODIS_ALG2(:,1)),1),'ro')
% hold off

% 
% aa = gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1);
% bb = GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1);
% 
% corr(gpp_day(~isnan(GPP_MODIS_ALG(:,1)),1),GPP_MODIS_ALG(~isnan(GPP_MODIS_ALG(:,1)),1))^2

% plot(170:1:299,gpp_day(:,1),'ro',170:1:299,GPP_MODIS_ALG(:,1),'bo')