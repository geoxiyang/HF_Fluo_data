%% Analyze diurnal patterns of GPP, SIF, LE, SIFy, LUE etc. in Harvard Forest 2013-2014
clc
clear variables


load('/Volumes/XiYangResearch/Projects/1.SCOPE_HF/4.matlab/hf_hourly_2013_newtest.mat',...
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
 
 GPP(GPP<0) = NaN;
 LE(LE<0)   = NaN;
 H(H<0)     = NaN;

 fit_scatter = cell(1,5);
 fit_mean    = cell(1,5);

 
 
 
 
 
 
 
% 
%  sif_hourly = sif_hourly';
% for month_i = 1:5
%     
%    lb = datenum(2013,month_i+5,1) - datenum(2013,1,1) +1;
%    ub = datenum(2013,month_i+6,1) - datenum(2013,1,1);
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
%  %  bin_width  = 0.1;
%    
%    % 30 bins
% %    for kk = 1:30
% %        subs = temp_sif >= (kk-1)*0.1 & temp_sif < kk*0.1 & temp_sp >0.5 & temp_sif > 0 & temp_GPP > 0 & temp_LE > 0;
% %        if sum(subs) == 0
% %           
% %            sif_bin_mean(kk,month_i)     = NaN;
% %            gpp_bin_mean(kk,month_i)     = NaN;
% %            le_bin_mean(kk,month_i)      = NaN;
% %            apar_bin_mean(kk,month_i)    = NaN;
% %            sif_bin_sd(kk,month_i)       = NaN;
% %            gpp_bin_sd(kk,month_i)       = NaN;
% %            le_bin_sd(kk,month_i)        = NaN;
% %            apar_bin_sd(kk,month_i)      = NaN;
% %            
% %        else
% % 
% %            sif_bin_mean(kk,month_i)     = nanmean(temp_sif(subs));
% %            gpp_bin_mean(kk,month_i)     = nanmean(temp_GPP(subs));
% %            le_bin_mean(kk,month_i)      = nanmean(temp_LE(subs));
% %            apar_bin_mean(kk,month_i)    = nanmean(temp_apar(subs));
% %            sif_bin_sd(kk,month_i)       = nanstd(temp_sif(subs));
% %            gpp_bin_sd(kk,month_i)       = nanstd(temp_GPP(subs));
% %            le_bin_sd(kk,month_i)        = nanstd(temp_LE(subs));
% %            apar_bin_sd(kk,month_i)      = nanstd(temp_apar(subs));  
% %            
% %        end  
% %    end
% 
%    for jj = 1:length(temp_doy1)-1
% 
%    % figure 1 all scatter 50% less cloud & temp_sp >0.5
%     positivesub = temp_sif>0 & temp_GPP>0  & temp_doy >= temp_doy1(jj) & temp_doy < temp_doy1(jj)+1;
%     xData       = temp_sif(positivesub);
%     yData       = temp_GPP(positivesub);
%     
%     if length(xData) <2 || length(yData) <2 
%         continue
%     end
%     
%     %figure
%     scatter(xData,yData,50,temp_apar(positivesub),'filled');
%     colorbar
%     hold on
%     
%     ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Robust = 'Off';
%     opts.StartPoint = [max(yData) 0.5];
%     opts.Lower  = [0 0];
%     opts.Upper  = [Inf Inf];
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     
%     h = plot(fitresult);
%     set(gca,'YLim',[0,50],...
%         'FontSize',16);
%     legend('off')
%     %legend( h, 'GPP vs. SIF', 'fit', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
%     ylabel('GPP(umol/m^{2}/sec)','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)
%     
% %    hh = text(0.1,max(yData)*1.1,['r^{2}=' num2str(gof.rsquare,'%4.2f')]);
%     hh.FontSize = 20;    
%          
%    end
% 
%     hold off
% 
%     %monthname = {'August';'September';'October'};
%     monthname = {'June';'July';'August';'September';'October'};
%     xl          = xlim;
%     yl          = ylim;
%     textpos     = [0.05*(xl(2)-xl(1)) + xl(1), yl(2) - 0.05*(yl(2)-yl(1))];
%     text(textpos(1),textpos(2),monthname{month_i});
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/monthly_SIF_GPP_dailyscatter_2013.ps';
%     
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end
%     
%     close(gcf);  


%    % figure 1 all scatter 50% less cloud
%     positivesub = temp_sif>0 & temp_LE>0 & temp_sp >0.5;
%     xData       = temp_sif(positivesub);
%     yData       = temp_LE(positivesub);
%     
%     figure
%     scatter(xData,yData,50,temp_apar(positivesub),'filled');
%     colorbar
%     hold on
%     
%     ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Robust = 'LAR';
%     opts.StartPoint = [max(yData) 0.5];
%     opts.Lower  = [0 0];
%     opts.Upper  = [2*max(yData) 5];
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     
%     h = plot(fitresult);
%     set(gca,'YLim',[0,max(yData)*1.2],...
%         'FontSize',16);
%     legend('off')
%     %legend( h, 'GPP vs. SIF', 'fit', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
%     ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)
%     
%     hh = text(0.1,max(yData)*1.1,['r^{2}=' num2str(gof.rsquare,'%4.2f')]);
%     hh.FontSize = 20;    
%     hold off
%     
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/monthly_SIF_LE_allscatter_2012.ps';
%     
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end
%     
%     close(gcf);        
    
    
%    % figure 2 mean scatter
%    
%     xData       = sif_bin_mean(~isnan(sif_bin_mean));
%     yData       = le_bin_mean(~isnan(sif_bin_mean));
%     xerr        = sif_bin_sd(~isnan(sif_bin_mean));
%     yerr        = le_bin_sd(~isnan(sif_bin_mean));
%     
%     figure
%     errorbarxy(xData,yData,xerr,yerr,'o');
%     hold on
%     
%     ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Robust = 'LAR';
%     opts.StartPoint = [max(yData) 0.5];
%     opts.Lower  = [0 0];
%     opts.Upper  = [2*max(yData) 5];
%     
%     % Fit model to data.
%     [fitresult2, gof] = fit( xData, yData, ft, opts );
%     
%     h = plot(fitresult2);
%     set(gca,'YLim',[0,max(yData)*1.2],...
%         'FontSize',16);
%     legend('off')
%     %legend( h, 'GPP vs. SIF', 'fit', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
%     ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)
%     
%     hh = text(0.1,max(yData)*1.1,['r^{2}=' num2str(gof.rsquare,'%4.2f')]);
%     hh.FontSize = 20;    
%     hold off
%     
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/monthly_SIF_LE_meanscatter_2013.ps';
%     
%     if exist(opfilename,'file')
%         print(gcf,'-dpsc',opfilename,'-append');
%     else
%         print(gcf,'-dpsc',opfilename);
%     end
%     
%     close(gcf); 
%     
%     % figure 3: compare the fits
%     
%     plot(fitresult);
%     hold on
%     plot(fitresult2);
%     
%     xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
%     ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)    
%     
%     opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/monthly_SIF_LE_fit_2013.ps';
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
%     fit_scatter{1,month_i} = fitresult;
%     fit_mean{1,month_i}    = fitresult2;
   
%    %only use sp > 0.5;
%    temp_sif(temp_sp<=0.5)       = NaN;
%    temp_apar(temp_sp<=0.5)      = NaN;
%    temp_pri2(temp_sp<=0.5)      = NaN;
%    temp_GPP(temp_sp<=0.5)       = NaN;
%    temp_LE(temp_sp<=0.5)        = NaN;   
%    
%    sif_mean_month_hourly(month_i,:)     = nanmean(temp_sif,2);
%    sif_sd_month_hourly(month_i,:)       = nanstd(temp_sif,0,2);
%    apar_mean_month_hourly(month_i,:)    = nanmean(temp_apar,2);
%    apar_sd_month_hourly(month_i,:)      = nanstd(temp_apar,0,2);
%    GPP_mean_month_hourly(month_i,:)     = nanmean(temp_GPP,2);
%    GPP_sd_month_hourly(month_i,:)       = nanstd(temp_GPP,0,2);
%    pri2_mean_month_hourly(month_i,:)    = nanmean(temp_pri2,2);
%    pri2_sd_month_hourly(month_i,:)      = nanstd(temp_pri2,0,2);
%    LE_mean_month_hourly(month_i,:)      = nanmean(temp_LE,2);
%    LE_sd_month_hourly(month_i,:)        = nanstd(temp_LE,0,2);
% 
% end

% colors = {'r','b','k','m','c'};
% for month_i = 1:5
%    
%     htemp = plot(fit_scatter{month_i});
%     set(htemp,'Color',colors{month_i},'LineWidth',2)
%     hold on
%     
% end
% 
% hold off
% 
% textleg = {'June','July','August','September','October'};
% legend(textleg)
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)    
% 
% opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/compare_SIF_LE_scatterfit_2013.png';
% print(gcf,'-dpng','-r300',opfilename);
% 
% close(gcf);
% 
% for month_i = 1:5
%    
%     htemp2 = plot(fit_mean{month_i});
%     set(htemp2,'Color',colors{month_i},'LineWidth',2)  
%     hold on
%     
%     
%     
% end
% 
% hold off
% legend(textleg)
% xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
% ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)    
% 
% opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/compare_SIF_LE_meanfit_2013.png';
% print(gcf,'-dpng','-r300',opfilename);
% close(gcf);

% sub = doy_hourly >=176 & doy_hourly <=177;
% 
% [ax,h1,h2]      = plotyy(doy_hourly(sub),GPP(sub),doy_hourly(sub),sif_hourly(sub));
% h1.LineStyle    = 'none';
% h2.LineStyle    = 'none';
% h1.Marker       = 'o';
% h2.Marker       = 'o';
% 
% sub1 =  doy_hourly >=176 & doy_hourly <=177 & GPP >0 & sif_hourly>0;
% figure
% scatter(sif_hourly(sub1),GPP(sub1));

sif_yield = sif_hourly./apar_hourly;
lue       = GPP./apar_hourly;

for ii = 5:5 %1:79 1:173 5:130
    
    sub1 =  doy_hourly >=ii+169 & doy_hourly <=ii+170 & GPP >0 & sif_hourly>0 & apar_hourly>0 & LE>0;
    
    %opfilepath = ['/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/diurnal/DOY_' num2str(ii+169) '.png'];

%    [xData1, yData1] = prepareCurveData(sif_hourly(sub1), GPP(sub1));
    [xData1, yData1] = prepareCurveData(sif_hourly(sub1), LE(sub1));

    [xData2, yData2] = prepareCurveData(lue(sub1), sif_yield(sub1));
    
    if length(xData1) < 6 || length(xData2) <6
        continue
    end
    
    
    VPD = vpd_hourly(sub1);
    
    % 1. PS file one: GPP-SIF, dots APAR colored
    
    figure
    
    scatter(xData1,yData1,50,apar_hourly(sub1),'filled');
%    scatter(xData1,yData1,50,sunportion_hourly(sub1),'filled');
    colorbar
    hold on
    
    ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Off';
    opts.StartPoint = [max(yData1) 0.5];
    opts.Lower  = [0 0];
    opts.Upper  = [2*max(yData1) 5];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData1, yData1, ft, opts );
    
    h = plot(fitresult);
    set(gca,'YLim',[0,max(yData1)*1.2],'XLim',[0,max(xData1)*1.2],...
        'FontSize',16);
    legend('off')
    %legend( h, 'GPP vs. SIF', 'fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel('SIF(mw/m^{2}/nm/sr)','FontSize',20);
    ylabel('LE(w/m^{2})','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)
    
    hh = text(0.1,max(yData1)*1.1,['DOY' num2str(ii+169) '   ' 'r^{2}=' num2str(gof.rsquare,'%4.2f')]);
    hh.FontSize = 20;    
    hold off
    
    opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/diurnal_LE_SIF_ROBUTSTOFF_APARcolored_summary_2012.ps';
    
    if exist(opfilename,'file')
        print(gcf,'-dpsc',opfilename,'-append');
    else
        print(gcf,'-dpsc',opfilename);
    end
    
    close(gcf);
    
    % 2. PS file two: GPP vs. APAR, SIF vs. APAR
    
    figure
    
    [ax,h1,h2] = plotyy(apar_hourly(sub1),GPP(sub1),apar_hourly(sub1),sif_hourly(sub1));
    h1.LineStyle = 'none';
    h2.LineStyle = 'none';
    h1.Marker    = '.';
    h2.Marker    = '.';
    h1.Color     = 'r';
    h2.Color     = 'b';
    h1.MarkerSize= 24;
    h2.MarkerSize= 24;
    
    set(ax(1),'XLim',[0,max(apar_hourly(sub1))*1.2],'YColor','r')
    set(ax(2),'XLim',[0,max(apar_hourly(sub1))*1.2],'YColor','b')
    
    hh = text(100,max(yData1)*0.9,['DOY' num2str(ii+126) '   ' 'r^{2}=' num2str(gof.rsquare,'%4.2f')]);
    hh.FontSize = 20; 

    set(ax,'FontSize',16);
    xlabel('APAR(umol/m^{2}/sec)','FontSize',20);
    ylabel(ax(1),'GPP(umol/m^{2}/sec)','FontSize',20);
    ylabel(ax(2),'SIF(mw/m^{2}/nm/sr)','FontSize',20);
    
    opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/diurnal_GPP_APAR_SIF_APAR_summary_2012.ps';
    
    if exist(opfilename,'file')
        print(gcf,'-dpsc',opfilename,'-append');
    else
        print(gcf,'-dpsc',opfilename);
    end
    
    close(gcf);    
    
    % 3. PS file three: SIFy vs. LUE
    
    figure
    scatter(xData2,yData2,50,apar_hourly(sub1),'filled');
    colorbar
    set(gca,'YLim',[0,max(yData2)*1.2],...
        'XLim',[0,max(xData2)*1.2],...
        'FontSize',16);
    hh = text(min(yData2)*1.1,max(yData2)*1.1,['DOY' num2str(ii+126)]);
    hh.FontSize = 20; 
    xlabel('LUE','FontSize',20);
    ylabel('SIFy','FontSize',20);%LE(w/m^{2})%GPP(umol/m^{2}/sec)    

    opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/diurnal_LUE_SIFy_APARcolored_summary_2012.ps';
    
    if exist(opfilename,'file')
        print(gcf,'-dpsc',opfilename,'-append');
    else
        print(gcf,'-dpsc',opfilename);
    end
    
    close(gcf); 
    
    
    % 4. PS file four: SIFy vs. APAR and LUE vs. APAR
    
    figure
    
    [ax,h1,h2] = plotyy(apar_hourly(sub1),xData2,apar_hourly(sub1),yData2);
    h1.LineStyle = 'none';
    h2.LineStyle = 'none';
    h1.Marker    = '.';
    h2.Marker    = '.';
    h1.Color     = 'r';
    h2.Color     = 'b';
    h1.MarkerSize= 24;
    h2.MarkerSize= 24;
    
    set(ax(1),'XLim',[0,max(apar_hourly(sub1))*1.2],'YColor','r')
    set(ax(2),'XLim',[0,max(apar_hourly(sub1))*1.2],'YColor','b')
    
    hh = text(100,max(yData1)*0.9,['DOY' num2str(ii+126) '   ' 'r^{2}=' num2str(gof.rsquare,'%4.2f')]);
    hh.FontSize = 20; 

    set(ax,'FontSize',16);
    xlabel('APAR(umol/m^{2}/sec)','FontSize',20);
    ylabel(ax(1),'LUE','FontSize',20);
    ylabel(ax(2),'SIFy','FontSize',20);
    
    opfilename = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/newfig/diurnal_LUE_APAR_SIFy_APAR_summary_2012.ps';
    
    if exist(opfilename,'file')
        print(gcf,'-dpsc',opfilename,'-append');
    else
        print(gcf,'-dpsc',opfilename);
    end
    
    close(gcf);  
    
    
    
    
    
% %    Set up fittype and options.
%     ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Robust = 'LAR';
%     opts.StartPoint = [max(yData) 0.5];
%     opts.Lower  = [0 0];
%     opts.Upper  = [2*max(yData) 5];
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     
%     figure( 'Name', 'untitled fit 1' );
%     h = plot( fitresult, xData, yData );
%     h(1).Marker     = '.';
%     h(1).MarkerSize = 20;
%     set(gca,'YLim',[0,max(yData)*1.2])
%     legend( h, 'GPP vs. SIF', 'fit', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel SIF(mw/m^{2}/nm/sr)
%     ylabel GPP(umol/m^{2}/sec)%LE(w/m^{2})%GPP(umol/m^{2}/sec)
%     
%     hh = text(0.1,max(yData)*1.1,['DOY' num2str(ii+169) '   ' 'r^{2}=' num2str(gof.rsquare,'%4.2f')]);
%     hh.FontSize = 20;
    

    
    
end







