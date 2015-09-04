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

work_dir    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_daily_2015-08-07-1150/';      % Change the input folder for the run you want to plot
ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/';
%time_dir    = '/Volumes/XiYangResearch/src/SCOPE/data/input/dataset HF_ts/';

wl          = dlmread([work_dir 'wl.dat'],'',2,0);
fEnergy     = dlmread([work_dir 'fluorescence.dat'],'',2,0);
fluxes      = dlmread([work_dir 'fluxes.dat'],'',2,0);
para_combination = dlmread([work_dir 'pars_and_input_short.dat'],'',1,0);

%test_sub    = para_combination(:,1) >= 40 & para_combination(:,1) <= 60 & para_combination(:,2) <= 80 &para_combination(:,2) >= 60 & para_combination(:,3) >= 4 & para_combination(:,3) <= 6 & para_combination(:,4) <= 1200;

LE_SCOPE    = fluxes(:,9);
A           = fluxes(:,11);
aPAR        = fluxes(:,18);
ftot        = fluxes(:,21);
ftotyield   = fluxes(:,22);
Rsp         = fluxes(:,16);
time        = fluxes(:,4);

GPPsc       = A + Rsp;

f740        = fEnergy(:,740-639);
% f737        = fEnergy(:,737-639);
f760        = fEnergy(:,760-639);                     % Common O2A retrieval wavelength
f755        = fEnergy(:,755-639);                     % at the GOSAT wavelenght. Unit is: W m-2 s-1 um-1 sr-1

timesub     = time>=175 & time<176;
ofile1      = [ofigure 'doy175_GPP_obs_mod.png'];
ofile2      = [ofigure 'doy175_SIF_obs_mod.png'];
ofile3      = [ofigure 'doy175_LUE_SIFy.png'];
ofile4      = [ofigure 'doy175_LUE_SIFy_diurnal.png'];
ofile5      = [ofigure 'seasonal_LUE_obs_mod.png'];

load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly.mat','GPP','H','LE','doy_hourly','sif_hourly','doy_ems','apar_hourly');

% plot(time,A./aPAR,'ro',doy_hourly,GPP./apar_hourly,'bo') %GPP./apar_hourly %A./aPAR
% ylim([0 0.1])
% print(gcf,'-dpng','-r300',ofile5);
% close(gcf)

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
% set(ax)

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
