%% SIF_Shadow.m
%  Work on the relationship between the shadow fraction calculated from
%  ArcMap Hillshade (Hilker et al., one of his RSE papers), and
%  SIF/SIFyield. We also compare the sunlit leaves fraction from the
%  estimation of FLiGHT model by Hideki Kobayashi. We focus on the data
%  from year 2013.

clear variables
clc

%% 1. Load all three data
% SIF data
load('/Volumes/XiYangBackUp/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newtimestamp.mat','sif_hourly','apar_hourly','sunportion_hourly');

% Shadow Fraction
shadowF = importdata('/Volumes/XiYangBackUp/Data/6.HFData/ShadowFraction/shadowfraction_latestrun.csv');
shadowF(shadowF<=0) = NaN;

% FLiGHT simulation
rawsim  = importdata('/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/2.Data/FLiGHT/SIFtot_summary2013.csv');
sunlai  = rawsim.data(:,3);
totlai  = rawsim.data(:,2);
doysim  = rawsim.data(:,1);

%% 2. Compare ShadowF with FLiGHT simulation

shalai  = totlai - sunlai;
shaF_mod= shalai./totlai;

shaF_mod= shaF_mod(doysim>=170&doysim<300,1);
shaF_hr = nan(3120,1);

for ii = 1:3120
    shaF_hr(ii,1) = nanmean(shaF_mod(((ii-1)*2+1):ii*2,1));
end

% subplot(2,1,1)
% plot(shadowF,shaF_hr,'ko')
% xlabel('LiDAR shadow fraction')
% ylabel('FLiGHT shadow fraction')
% 
% subplot(2,1,2)
% plot(1:3120,shadowF,'k^',1:3120,shaF_hr,'rx')
% 
% print('/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/shadowF_comp.eps',gcf,'-depsc','-tiff')
% close(gcf)

%% 3. Compare SIFyield with Shadow Fractions

SIFyield_hourly = sif_hourly./apar_hourly;
SIFyield_hourly(SIFyield_hourly==0) = NaN;

% figure
% plot(shadowF,SIFyield_hourly,'ko')
% xlabel('Shadow Fraction From LiDAR')
% ylabel('SIFyield')
% print('/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/ShadowFractionLiDAR_vs_SIFyield.eps',gcf,'-depsc','-tiff')
% close(gcf)
% 
% figure
% plot(shadowF,SIFyield_hourly,'ko')
% xlabel('Shadow Fraction From LiDAR')
% ylabel('SIFyield')
% ylim([0 0.005])
% print('/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/ShadowFractionLiDAR_vs_SIFyield_ENLARGED.eps',gcf,'-depsc','-tiff')
% close(gcf)
% 
% 
% figure
% plot(shaF_hr,SIFyield_hourly,'ro')
% xlabel('Shadow Fraction From FLiGHT')
% ylabel('SIFyield')
% print('/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/ShadowFractionFLiGHT_vs_SIFyield.eps',gcf,'-depsc','-tiff')
% close(gcf)

for jj = 1:130
    
    shadowF_daily(jj,1) = nanmean(shadowF(((jj-1)*24+1):jj*24,1));
    shaF_daily(jj,1)    = nanmean(shaF_hr(((jj-1)*24+1):jj*24,1));
    SIFyield_daily(jj,1)= nanmean(SIFyield_hourly(((jj-1)*24+1):jj*24,1));
    sp_daily(jj,1)      = nanmean(sunportion_hourly(((jj-1)*24+1):jj*24,1));
    
    
    
%     figure
%     scatter(shadowF(((jj-1)*24+1):jj*24,1),SIFyield_hourly(((jj-1)*24+1):jj*24,1),64,sunportion_hourly(((jj-1)*24+1):jj*24,1))
%     title(num2str(jj))
%     colorbar
%     
%     opfilename = '/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/shadowF_LiDAR_SIFyield_daily.ps';
% 
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end  
%     close(gcf)
%     
%     figure
%     scatter(shaF_hr(((jj-1)*24+1):jj*24,1),SIFyield_hourly(((jj-1)*24+1):jj*24,1),64,sunportion_hourly(((jj-1)*24+1):jj*24,1))
%     title(num2str(jj))
%     colorbar
%     
%     opfilename = '/Volumes/XiYangBackUp/Projects/10.HF_SIF_Synthesis/1.JPG/shadowF_FLiGHT_SIFyield_daily.ps';
% 
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end      
%     close(gcf)
    
    
end

shadowF_daily(sp_daily<=0.2) = NaN;
SIFyield_daily(sp_daily<=0.2) = NaN;

yyaxis left
plot(1:130,shadowF_daily,'ko')
yyaxis right
plot(1:130,SIFyield_daily,'rx')





