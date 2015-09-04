%% more figures to make

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
             'ea_hourly',...
             'rh_hourly');
         
data2013 = load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newcutoff.mat',...
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
             'ea_hourly',...
             'rh_hourly');

data2012.vpd_hourly = data2012.ea_hourly .* (100-data2012.rh_hourly)./data2012.rh_hourly;
data2014.vpd_hourly = data2014.ea_hourly .* (100-data2014.rh_hourly)./data2014.rh_hourly;

         
         
gpp_hourly_all   = [data2012.GPP;data2013.GPP;data2014.GPP];
le_hourly_all    = [data2012.LE;data2013.LE;data2014.LE];
apar_hourly_all  = [data2012.apar_hourly;data2013.apar_hourly;data2014.apar_hourly];
sif_hourly_all   = [data2012.sif_hourly;data2013.sif_hourly;data2014.sif_hourly];
doy_hourly_all   = [data2012.doy_hourly;data2013.doy_hourly+366;data2014.doy_hourly+366+365];
sp_hourly_all    = [data2012.sunportion_hourly;data2013.sunportion_hourly;data2014.sunportion_hourly];

th = 0.0;

gpp_hourly_all(gpp_hourly_all<=0 | sp_hourly_all<th)     = NaN;
le_hourly_all(le_hourly_all<=0 | sp_hourly_all<th)       = NaN;
apar_hourly_all(apar_hourly_all<=0 | sp_hourly_all<th)   = NaN;
sif_hourly_all(sif_hourly_all<=0 | sp_hourly_all<th)     = NaN;

%% Hourly relationship between GPP and SIF; GPP and LE.

% dataset = data2014;
% 
% for ii=1:11
%    
%     time = dataset.doy_hourly - floor(dataset.doy_hourly);
%     
%     temp_sub = time >= (ii+4.9)/24.0 & time <= (ii+5.1)/24.0 & dataset.sif_hourly>0 & dataset.sunportion_hourly>th;
%     temp_gpp = dataset.GPP(temp_sub);
%     temp_sif = dataset.sif_hourly(temp_sub);
%     temp_le  = dataset.LE(temp_sub);
%     temp_apar= dataset.apar_hourly(temp_sub);
%     
%     scatter(temp_sif,temp_gpp,50,temp_apar,'filled');
%     colorbar
%     xl = xlim;
%     yl = ylim;
%     
%     
%     text(0.1*(xl(2)-xl(1))+xl(1),0.9*(yl(2)-yl(1))+yl(1),num2str(ii+5),'FontSize',24);
%    
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/fixedhour_SIF_GPP_SP0_2014.ps';
%     
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end
%     
%     close(gcf);
%     
%     %====================
%     scatter(temp_sif,temp_le,50,temp_apar,'filled');
%     colorbar
%     xl = xlim;
%     yl = ylim;
%     
%     
%     text(0.1*xl(2),0.9*yl(2),num2str(ii+5),'FontSize',24);
%    
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/fixedhour_SIF_LE_SP0_2014.ps';
%     
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end
%     
%     close(gcf);
%        
%     
% end

%% Seasonal relationship summary, and check if it is robust

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
vpd_daily_2012  = nan(length(doy_tmp_2012),1);
vpd_daily_2013  = nan(length(doy_tmp_2013),1);
vpd_daily_2014  = nan(length(doy_tmp_2014),1);

for ii = 1:length(doy_tmp_2012)
    
    gpp_daily_2012(ii,1)    = nanmean(data2012.GPP(data2012.GPP >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<(doy_tmp_2012(ii)+1)));
    le_daily_2012(ii,1)     = nanmean(data2012.LE(data2012.LE >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    sif_daily_2012(ii,1)    = nanmean(data2012.sif_hourly(data2012.sif_hourly >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
    apar_daily_2012(ii,1)   = nanmean(data2012.apar_hourly(data2012.apar_hourly>0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));
%    gpp_max_2012(ii,1)      = nanmax(data2012.GPP(data2012.GPP >0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<(doy_tmp_2012(ii)+1)));
    vpd_daily_2012(ii,1)   = nanmean(data2012.vpd_hourly(data2012.vpd_hourly>0 & data2012.doy_hourly>=doy_tmp_2012(ii) & data2012.doy_hourly<doy_tmp_2012(ii)+1));

end

for ii = 1:length(doy_tmp_2013)
    
    gpp_daily_2013(ii,1)    = nanmean(data2013.GPP(data2013.GPP >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    le_daily_2013(ii,1)     = nanmean(data2013.LE(data2013.LE >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    sif_daily_2013(ii,1)    = nanmean(data2013.sif_hourly(data2013.sif_hourly >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
    apar_daily_2013(ii,1)   = nanmean(data2013.apar_hourly(data2013.apar_hourly>0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));
%    gpp_max_2013(ii,1)      = nanmax(data2013.GPP(data2013.GPP >0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<(doy_tmp_2013(ii)+1)));
    vpd_daily_2013(ii,1)   = nanmean(data2013.vpd_hourly(data2013.vpd_hourly>0 & data2013.doy_hourly>=doy_tmp_2013(ii) & data2013.doy_hourly<doy_tmp_2013(ii)+1));

end

for ii = 1:length(doy_tmp_2014)
    
    gpp_daily_2014(ii,1)    = nanmean(data2014.GPP(data2014.GPP >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    le_daily_2014(ii,1)     = nanmean(data2014.LE(data2014.LE >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    sif_daily_2014(ii,1)    = nanmean(data2014.sif_hourly(data2014.sif_hourly >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
    apar_daily_2014(ii,1)   = nanmean(data2014.apar_hourly(data2014.apar_hourly>0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));
%    gpp_max_2014(ii,1)      = nanmax(data2014.GPP(data2014.GPP >0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<(doy_tmp_2014(ii)+1)));
    vpd_daily_2014(ii,1)   = nanmean(data2014.vpd_hourly(data2014.vpd_hourly>0 & data2014.doy_hourly>=doy_tmp_2014(ii) & data2014.doy_hourly<doy_tmp_2014(ii)+1));

end




gpp_daily_all   = [gpp_daily_2012;gpp_daily_2013;gpp_daily_2014];
le_daily_all    = [le_daily_2012;le_daily_2013;le_daily_2014];
apar_daily_all  = [apar_daily_2012;apar_daily_2013;apar_daily_2014];
sif_daily_all   = [sif_daily_2012;sif_daily_2013;sif_daily_2014];
doy_daily_all   = [doy_tmp_2012;doy_tmp_2013+366;doy_tmp_2014+366+365];
gpp_max_all     = [gpp_max_2012;gpp_max_2013;gpp_max_2014];


%plot(doy_daily_all,sif_daily_all,'ro')

%gpp_daily_2013(gpp_daily_2013>20) = NaN;
sif_daily_2014(sif_daily_2014>1) = NaN;

plot(sif_daily_2012,gpp_daily_2012,'ro');
hold on

plot(sif_daily_2013,gpp_daily_2013,'bo');
plot(sif_daily_2014,gpp_daily_2014,'ko');

%lsline


ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Off';
opts.StartPoint = [max(gpp_daily_all) 0.5];
opts.Lower  = [0 0];
opts.Upper  = [2*max(gpp_daily_all) 5];

% Fit model to data.
[fitresult, gof]     = fit( sif_daily_2012(~isnan(sif_daily_2012) & ~isnan(gpp_daily_2012)), gpp_daily_2012(~isnan(sif_daily_2012) & ~isnan(gpp_daily_2012)), ft, opts );
[fitresult2, gof2]   = fit( sif_daily_2013(~isnan(sif_daily_2013) & ~isnan(gpp_daily_2013)), gpp_daily_2013(~isnan(sif_daily_2013) & ~isnan(gpp_daily_2013)), ft, opts );
[fitresult3, gof3]   = fit( sif_daily_2014(~isnan(sif_daily_2014) & ~isnan(gpp_daily_2014)), gpp_daily_2014(~isnan(sif_daily_2014) & ~isnan(gpp_daily_2014)), ft, opts );

h   = plot(fitresult);
h2  = plot(fitresult2);
h3  = plot(fitresult3);
set(h,'Color','r','LineWidth',2)
set(h2,'Color','b','LineWidth',2)
set(h3,'Color','k','LineWidth',2)


set(gca,'YLim',[0,max(gpp_daily_all)*1.5],...
    'FontSize',16);
legend({'2012';'2013';'2014'})
xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
ylabel('GPP(umol/m^{2}/sec)','FontSize',20);

hold off

print(gcf,'-dpng','-r300','/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/sif_gpp_diffseasons_sif2014LG1.png')
close(gcf)




