%% fluo_analysis

% Author: Xi Yang (geoxiyang@gmail.com)

% Date: Aug. 8, 2014

% Read the SCOPE-produced fluorescence file
% Plot fluorescence with one parameter changes and the rest fixed
% Calculate the ratio between fluorescence at a single wavelength and 
% the integral of SIF

% In the future, this file make a LUT, and use it as the input to SCOPE,
% and then we plot the response of fluorescence as a function of variables

clear all
clc

global constants
[constants] = define_constants();
 %HF_ts_daily_2015-09-15-1116
work_dir    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_daily_2015-10-07-1544/';      % Change the input folder for the run you want to plot
ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/';
%time_dir    = '/Volumes/XiYangResearch/src/SCOPE/data/input/dataset HF_ts/';

load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newtest.mat','GPP','H','LE','doy_hourly','sif_hourly','doy_ems','apar_hourly','sunportion_hourly');

wl          = dlmread([work_dir 'wl.dat'],'',2,0);
fEnergy     = dlmread([work_dir 'fluorescence.dat'],'',2,0);
fluxes      = dlmread([work_dir 'fluxes.dat'],'',2,0);
para_combination = dlmread([work_dir 'pars_and_input_short.dat'],'',1,0);
ref         = dlmread([work_dir 'reflectance.dat'],'',2,0);

ndvi_SCOPE  = (ref(:,750-400+1)-ref(:,685-400+1))./(ref(:,750-400+1) + ref(:,685-400+1));
pri_SCOPE   = (ref(:,531-400+1)-ref(:,570-400+1))./(ref(:,531-400+1) + ref(:,570-400+1));
%test_sub    = para_combination(:,1) >= 40 & para_combination(:,1) <= 60 & para_combination(:,2) <= 80 &para_combination(:,2) >= 60 & para_combination(:,3) >= 4 & para_combination(:,3) <= 6 & para_combination(:,4) <= 1200;

LE_SCOPE    = fluxes(:,9);
A           = fluxes(:,11);
aPAR        = fluxes(:,18);
ftot        = fluxes(:,21);
ftotyield   = fluxes(:,22);
Rsp         = fluxes(:,16);
time        = fluxes(:,4);

f740        = fEnergy(:,740-639);
% f737        = fEnergy(:,737-639);
f760        = fEnergy(:,760-639);                     % Common O2A retrieval wavelength
f755        = fEnergy(:,755-639);                     % at the GOSAT wavelenght. Unit is: W m-2 s-1 um-1 sr-1

% doy_SCOPE   = unique(fix(time));
% doy_obs     = unique(fix(doy_hourly));
% 
% for ii = 1:length(doy_SCOPE)
%     
%     lb = doy_SCOPE(ii) + 0.25;
%     ub = doy_SCOPE(ii) + 0.75;
%     daily_model_A(ii)       = nanmean(A(time>=lb & time<ub));
%     daily_model_LE(ii)      = nanmean(LE_SCOPE(time>=lb & time<ub));
%     
% end
% 
% for ii = 1:length(doy_SCOPE)
%     
%     lb = doy_SCOPE(ii) + 0.42;
%     ub = doy_SCOPE(ii) + 0.60;
%     daily_model_sif(ii)     = nanmean(f760(time>=lb & time<ub));
%     
% end
% 
% for jj = 1:length(doy_obs)
%    
%     lb = doy_obs(jj) + 0.25;
%     ub = doy_obs(jj) + 0.75;
%     daily_obs_GPP(jj) = nanmean(GPP(doy_hourly>=lb & doy_hourly<ub));
%     daily_obs_LE(jj)  = nanmean(LE(doy_hourly>=lb & doy_hourly<ub));
% 
% end
% 
% for jj = 1:length(doy_obs)
%    
%     lb = doy_obs(jj) + 0.42;
%     ub = doy_obs(jj) + 0.60;
%     daily_obs_sif(jj) = nanmean(sif_hourly(doy_hourly>=lb & doy_hourly<ub));
% 
% end
% 
% doy_obs_1 = doy_obs(6:127);
% daily_obs_sif = daily_obs_sif(6:127);
% daily_obs_GPP = daily_obs_GPP(6:127);
% daily_obs_LE  = daily_obs_LE(6:127);
% sub_sif = daily_obs_sif >0;
% 
% 
% subplot(3,1,1)
% plot(doy_SCOPE,daily_model_A,'ro',doy_obs_1,daily_obs_GPP,'bo')
% title('GPP')
% subplot(3,1,2)
% plot(doy_SCOPE,daily_model_sif,'ro',doy_obs_1(sub_sif),daily_obs_sif(sub_sif),'bo')
% title('SIF')
% subplot(3,1,3)
% plot(doy_SCOPE,daily_model_LE,'ro',doy_obs_1,daily_obs_LE,'bo')
% title('LE')
% print(gcf,'-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/model_data_comp.png')
% close(gcf)


day         = 175;
timesub     = time>=day & time<day+1;
ofile       = [ofigure 'doy' num2str(day) '_obs_mod_fluspect_2014_lowfqe0015_tto30.png'];
ofile1      = [ofigure 'doy' num2str(day) '_GPP_obs_mod.png'];
ofile2      = [ofigure 'doy' num2str(day) '_SIF_obs_mod.png'];
ofile3      = [ofigure 'doy' num2str(day) '_LE_obs_mod.png'];
% ofile3      = [ofigure 'doy175_LUE_SIFy.png'];
% ofile4      = [ofigure 'doy175_LUE_SIFy_diurnal.png'];
% ofile5      = [ofigure 'seasonal_LUE_obs_mod.png'];

time1.year  = 2013;
time1.month  = 6;
time1.day    = 23;
time1.min    = 0;
time1.sec    = 0;
time1.UTC    = -5;
location.latitude = 42.5;
location.longitude = -72.2;
location.altitude = 100; 


for ii = 1:24
    
    time1.hour = ii;
    
    sun_pos = sun_position(time1,location);
    
    sun_azimuth(ii) = sun_pos.azimuth;
    sun_zenith(ii)  = sun_pos.zenith;
    
end

% plot(1:1:24,sun_azimuth,'k*')
% hold on
% plot(1:1:24,sun_zenith,'ko')
% hold off
% 
% plot(time(timesub),pri_SCOPE(timesub),'bo')
% xlabel('hours')
% ylabel('pri')

%cos(sun_zenith'.*pi()/180.0)     %./cos(60.*pi()/180.0)


sub0  = doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub));


figure
plot(time(timesub),f760(timesub),'r*',doy_ems(sub0),sif_hourly(sub0),'k*');
figure
plot(time(timesub),f760(timesub),'r*',doy_ems(sub0),sif_hourly(sub0)./cos(sun_zenith'.*pi()/180.0),'k*');

figure
plot(time(timesub),f760(timesub)./aPAR(timesub),'ro-',doy_ems(sub0),sif_hourly(sub0)./(apar_hourly(sub0).*cos(sun_zenith'.*pi()/180.0)),'ko-');

figure
plot(time(timesub),A(timesub)./aPAR(timesub),'ro-',doy_ems(sub0),GPP(sub0)./(apar_hourly(sub0)),'ko-')

figure
plot(GPP(sub0)./(apar_hourly(sub0).*cos(sun_zenith'.*pi()/180.0)),sif_hourly(sub0)./(apar_hourly(sub0).*cos(sun_zenith'.*pi()/180.0)),'mo')
hold on
plot(A(timesub)./aPAR(timesub),f760(timesub)./aPAR(timesub),'b*')
hold off

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])

subplot(1,3,1)
plot(time(timesub),A(timesub),'ro',doy_ems(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),GPP(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),'ko');
title('GPP','FontSize',16)
% print(gcf,'-dpng','-r300',ofile5);
% close(gcf)
% 
subplot(1,3,2)
plot(time(timesub),f760(timesub),'r*',doy_ems(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),sif_hourly(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),'k*') %GPP./apar_hourly %A./aPAR
title('SIF','FontSize',16)

% ylim([0 0.1])
% print(gcf,'-dpng','-r300',ofile5);
% close(gcf)
% 
subplot(1,3,3)
plot(time(timesub),LE_SCOPE(timesub),'r^',doy_ems(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),LE(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),'k^') %GPP./apar_hourly %A./aPAR
title('LE','FontSize',16)

set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,6])
print(gcf,'-dpng','-r300',ofile);
close(gcf)

% ylim([0 0.1])
% print(gcf,'-dpng','-r300',ofile5);
% close(gcf)
% 
% figure
% plot(time(timesub),LE_SCOPE(timesub),'bo',doy_ems(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),LE(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),'ko');



% [ax,h1,h2] = plotyy(time,A,doy_hourly,GPP);
% h1.LineStyle = 'none';
% h2.LineStyle = 'none';
% h1.Marker    = '.';
% h2.Marker    = '.';
% h1.Color     = 'r';
% h2.Color     = 'b';
% h1.MarkerSize= 24;
% h2.MarkerSize= 24;
% set(ax(1),'YColor','r')

% figure
% plot(time(timesub),A(timesub),'ro',doy_ems(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),GPP(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))),'ko');
% print(gcf,'-dpng','-r300',ofile1);
% close(gcf)
% 
% 
% figure
% plot(time(timesub),f760(timesub),'rx',doy_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub))),sif_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub))),'kx')
% print(gcf,'-dpng','-r300',ofile2);
% close(gcf)
% 
% figure
% 
% timesub1 =  A(timesub)>0 & f760(timesub)>0 & GPP(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub))) >0 & sif_hourly(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub)))>0;
% time1    = time(timesub);
% time2    = doy_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub)));
% LUE_mod  = A(timesub)./aPAR(timesub);
% SIFy_mod = f760(timesub)./aPAR(timesub);
% LUE_obs  = GPP(doy_ems>=min(time(timesub)) & doy_ems<=max(time(timesub)))./apar_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub)));
% SIFy_obs = sif_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub)))./apar_hourly(doy_hourly>=min(time(timesub)) & doy_hourly<=max(time(timesub)));
% 
% subplot(2,1,1)
% scatter(LUE_mod(timesub1),SIFy_mod(timesub1),40,'filled')
% xlim([0 0.05])
% 
% subplot(2,1,2)
% scatter(LUE_obs,SIFy_obs,40,'filled')
% xlim([0 0.05])
% print(gcf,'-dpng','-r300',ofile3);
% close(gcf)
% 
% figure
% subplot(2,1,1)
% plotyy(time1(timesub1),LUE_mod(timesub1),time1(timesub1),SIFy_mod(timesub1))
% subplot(2,1,2)
% plotyy(time2(timesub1),LUE_obs(timesub1),time2(timesub1),SIFy_obs(timesub1))
% print(gcf,'-dpng','-r300',ofile4);
% close(gcf)
% % plot(time(1:4120,1),f760,'bo')
% sif_day = zeros(254-170+1,1);
% 
% for ii = 1:254-170
%     lb = ii-1+170;
%     ub = ii+170;
%     
%     time_sub = time>=lb & time <ub;
%     sif_day(ii,1)  = nanmean(f760(time_sub));
%     
%     
%     
%     
% end
% plot(170:254,sif_day,'bo')





% for ii = 1:size(f740)
% %     fPhotons    = 1E6   * e2phot(1E-9*wl(wl>639 & wl<851),pi*fEnergy(ii,:));   % convert into umol m-2 s-1 um-1
% %     fPhotonsInt = 0.001 * Sint(fPhotons,wl(wl>639 & wl<851));    % convert into umol m-2 s-1
% 
%     Factor(ii)      = sum(fEnergy(ii,:))/f740(ii);
%     Factor2(ii)      = sum(fEnergy(ii,:))/f755(ii);
% %     Factor2(ii)     = f760(ii)/f737(ii);
% %     Factor3(ii)     = f755(ii)/f737(ii);
% %     Factor4(ii)     = f760(ii)/f740(ii);
% end
% 
% 
% para_combination = dlmread([work_dir 'pars_and_input_short.dat'],'',1,0);
% plot(para_combination,Factor','r*-');
% 
% xlabel('V_{cmax}(umol m^{-2} s^{-1})','FontSize',16,'FontName','Helvetica');
% ylabel('SIF/F740','FontSize',16,'FontName','Helvetica');
% set(gca,'FontSize',16,'FontName','Helvetica');
% print('-dpng','-r300',[ofigure 'factor-vcmax-SIF740factor.png']);
% %print('-depsc','-r300',[ofigure 'factor-vcmax-chl40.eps']);


%% Linear regression
% Between Factors and independent varibles vcmax etc.
% x = [ones(size(para_combination(:,1))),para_combination];
% 
% [b,bint,r,rint,stats] = regress(Factor',x);
% 
% b
% stats
% 
% plot(b'*x',Factor,'bo',min(Factor):(max(Factor)-min(Factor))/100.0:max(Factor),min(Factor):(max(Factor)-min(Factor))/100.0:max(Factor),'r-')
% xlabel('Regression F755/SIF','FontSize',16,'FontName','Helvetica');
% ylabel('SCOPE F755/SIF','FontSize',16,'FontName','Helvetica');
% set(gca,'FontSize',16,'FontName','Helvetica');
% print('-dpng','-r300',[ofigure 'factor-vcmax-chl-regress_newSCOPE.png']);
% print('-depsc','-r300',[ofigure 'factor-vcmax-chl-regress_newSCOPE.eps']);



%% Draw figures
% rawdata = [Factor',para_combination];
% [pi,pj] = size(rawdata);
% 
% uni_para = unique(rawdata(:,2));
% fig = figure;
% cc  = hsv(10);
% M = zeros(size(uni_para));
% 
% xlabel('V_{cmax}(umol m^{-2} s^{-1})','FontSize',16,'FontName','Helvetica');
% ylabel('F755/SIF','FontSize',16,'FontName','Helvetica');
% xlim([0,130]);
% set(gca,'FontSize',16,'FontName','Helvetica');
% 
% for jj = 1:size(uni_para)
%     ind  = find(rawdata(:,2) == uni_para(jj));
%     data2 = rawdata(ind,:);
%     hold on
%     plot(data2(:,3),data2(:,1),'-*','Color',cc(jj,:));
%     M(jj) = uni_para(jj);
% end
% hold off
% outM = strtrim(cellstr(num2str(M))');
% 
% legend(gca,outM{:});
% print(gcf,'-dpng','-r300',[ofigure 'factor-vcmax-chl_newSCOPE.png']);
% print(gcf,'-depsc','-r300',[ofigure 'factor-vcmax-chl_newSCOPE.eps']);                                          
