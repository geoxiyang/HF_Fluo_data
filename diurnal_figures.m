% Making figures of diurnal patterns of
% SIF
% GPP


 load('SIF760daily.mat','halfhourly_result','raw_final_result');
 load('HF_2013_GPP.mat'); 
%  load('SIF760_result.mat');
 load('hf_barn_2013_env.mat','apar','fapar','par');

%  mkdir('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\','sif_gpp\');
%  mkdir('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\','sif_gpp_linear\');
%  mkdir('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\','sifyield\');
%  mkdir('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\','sif_gpp_mm\');
%  mkdir('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\','sifyield_linear\');
 
 %SIF-GPP linear
rsq_total = nan(numel(gpp_doy),1);
p_total = nan(numel(gpp_doy),1);
rmse_total = nan(numel(gpp_doy),1);
%SIF-GPP Michaelis-Menten
rsq1_total = nan(numel(gpp_doy),1);
p1_total = nan(numel(gpp_doy),1);
rmse1_total = nan(numel(gpp_doy),1);
%SIFyield-GPP linear
rsq2_total = nan(numel(gpp_doy),1);
p2_total = nan(numel(gpp_doy),1);
rmse2_total = nan(numel(gpp_doy),1);

 for doy_i=81:81 %1:numel(gpp_doy)


   lb = double(gpp_doy(doy_i))+0.;
   ub = double(gpp_doy(doy_i))+1.;
   temp_time = 24.0*(doy(doy>= lb & doy <ub) - gpp_doy(doy_i));
   temp_sif = halfhourly_result(halfhourly_result(:,1)>= lb & halfhourly_result(:,1)<ub,2);
   temp_gpp = gpp_raw(doy>= lb & doy <ub);
   temp_apar = apar(doy>= lb & doy <ub);
   temp_par = par(doy>= lb & doy <ub);
   temp_yield = temp_sif./temp_apar;
   temp_lue = temp_gpp./temp_apar;
   
   temp_time_1 = 24.0*(raw_final_result(raw_final_result(:,1) >= lb & raw_final_result(:,1) < ub & raw_final_result(:,3) >=0.99, 1) - gpp_doy(doy_i));
   temp_array_1 = raw_final_result(raw_final_result(:,1) >= lb & raw_final_result(:,1) < ub & raw_final_result(:,3) >=0.99, 2);
   
   gpp_pos_sub = temp_gpp>0;
   temp_gpp = temp_gpp(gpp_pos_sub);
   temp_sif = temp_sif(gpp_pos_sub);
   temp_time= temp_time(gpp_pos_sub);
   temp_yield = temp_yield(gpp_pos_sub);
   temp_apar = temp_apar(gpp_pos_sub);
   
   
   %Figure 1   
   fig = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig,'Visible','off');
   [AX,H1,H2] = plotyy(temp_time,temp_sif,temp_time,temp_gpp);
   set(H1,'marker','O');
   set(H2,'marker','O');
   set(AX,'xlim',[0,24]);
   set(AX(1),'box','off'); % Primary Y axis ticks only show on one side
   set(get(AX(1),'Ylabel'),'String','SIF(mw/m2/sr/nm)','FontSize',15,'FontName','Whitney-Book');
   set(get(AX(2),'Ylabel'),'String','GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   set(AX(1),'YLim',[0 3],'FontSize',12,'FontName','Whitney-Book');
   set(AX(2),'YLim',[0 50],'FontSize',12,'FontName','Whitney-Book');
   
   xlabel('Hours','FontSize',16,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   text(pos(1),pos(2),num2str(gpp_doy(doy_i)),'FontSize',20,'FontName','Whitney-Book');
   %=================
   
   %Figure 2
   pp = polyfit(temp_sif(temp_sif>0),temp_gpp(temp_sif>0),1);
   fitresult = linearfit(polyval(pp,temp_sif(temp_sif>0)),temp_gpp(temp_sif>0));
   rsq = fitresult(6);
   pvalue = fitresult(7);
   slope = fitresult(1);
   intercept = fitresult(3); 
   rmse = fitresult(5);
   
   rsq_total(doy_i) = rsq;
   p_total(doy_i) = pvalue;
   rmse_total(doy_i) = rmse;   
   
   fig2 = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig2,'Visible','off');
   start = min(temp_sif(temp_sif>0));
   stop = max(temp_sif(temp_sif>0))+0.1*max(temp_sif(temp_sif>0)); %pad 10%
   plot(temp_sif(temp_sif>0),temp_gpp(temp_sif>0),'ro',...
       start:0.1:stop,polyval(pp,start:0.1:stop),'r-');
   axis([0 3 0 50]);
   
   xlabel('SIF(mw/sr/m2/nm)','FontSize',15,'FontName','Whitney-Book');   
   ylabel('GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   info = [num2str(gpp_doy(doy_i),'%d'),'  ','R2=',num2str(rsq,' %4.2f'),' ','RMSE=',num2str(rmse,' %4.2f')];
   text(pos(1),pos(2),info,'FontSize',20,'FontName','Whitney-Book');
   pos2 = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.15*ht];
   eq = ['GPP=' num2str(pp(2),'%6.4f') '+' num2str(pp(1),'%6.4f') '*SIF   ' 'p=' num2str(pvalue,'%6.4f') ];
   text(pos2(1),pos2(2),eq ,'FontSize',15,'FontName','Whitney-Book');
   %=================

   %Figure 3
   
   fig3 = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig3,'Visible','Off');
   [AX,H1,H2] = plotyy(temp_time,temp_yield,temp_time,temp_gpp);
   set(H1,'marker','O');
   set(H2,'marker','O');
   set(AX,'xlim',[0,24]);
   set(AX(1),'box','off'); % Primary Y axis ticks only show on one side
   set(get(AX(1),'Ylabel'),'String','SIFyield','FontSize',15,'FontName','Whitney-Book');
   set(get(AX(2),'Ylabel'),'String','GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   set(AX(1),'YLim',[0 0.001],'FontSize',12,'FontName','Whitney-Book');
   set(AX(2),'YLim',[0 50],'FontSize',12,'FontName','Whitney-Book');

   xlabel('Hours','FontSize',16,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   text(pos(1),pos(2),num2str(gpp_doy(doy_i)),'FontSize',20,'FontName','Whitney-Book');    
   
   %=================
   
   %Figure 4
   x0 = [20 1];
   lb_mm = [-inf,-inf];
   ub_mm = [inf,inf];
   options = optimset('TolX',1e-10,'TolFun',1e-10,'Display','off');
   [x,SSresid] = lsqcurvefit(@michaelis_menten,x0,temp_sif(temp_sif>0),temp_gpp(temp_sif>0),lb_mm,ub_mm,options);
    
   SStotal = (length(temp_gpp(temp_sif>0))-1) * var(temp_gpp(temp_sif>0));
   rsq1_total(doy_i) = 1 - SSresid/SStotal;
   rmse1_total(doy_i) = (SSresid/(length(temp_gpp(temp_sif>0))))^(0.5);
   
   fig4 = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig4,'Visible','off');
   start = min(temp_sif(temp_sif>0));
   stop = max(temp_sif(temp_sif>0))+0.1*max(temp_sif(temp_sif>0)); %pad 10%
   plot(temp_sif(temp_sif>0),temp_gpp(temp_sif>0),'ro',...
       start:0.1:stop,michaelis_menten(x,transpose(start:0.1:stop)),'r-');
   axis([0 3 0 50]);
   
   xlabel('SIF(mw/sr/m2/nm)','FontSize',15,'FontName','Whitney-Book');   
   ylabel('GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   info = [num2str(gpp_doy(doy_i),'%d'),'  ','R2=',num2str(rsq1_total(doy_i),' %4.2f'),' ','RMSE=',num2str(rmse1_total(doy_i),' %4.2f')];
   text(pos(1),pos(2),info,'FontSize',20,'FontName','Whitney-Book');
   pos2 = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.15*ht];
   eq = ['GPP=' num2str(x(1),'%6.4f') '*SIF/(' num2str(x(2),'%6.4f') '+SIF)'];
   text(pos2(1),pos2(2),eq ,'FontSize',15,'FontName','Whitney-Book');
   
   %==========================
   
     %Figure 5
   pp5 = polyfit(temp_yield(temp_sif>0),temp_gpp(temp_sif>0),1);
   fitresult = linearfit(polyval(pp5,temp_yield(temp_sif>0)),temp_gpp(temp_sif>0));
   rsq2 = fitresult(6);
   pvalue2 = fitresult(7);
   slope2 = fitresult(1);
   intercept2 = fitresult(3);   
   rmse2 = fitresult(5);

   rsq2_total(doy_i) = rsq2;
   p2_total(doy_i) = pvalue2;
   rmse2_total(doy_i) = rmse2;   
   
   fig5 = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig5,'Visible','off');
   start = min(temp_yield(temp_sif>0));
   stop = max(temp_yield(temp_sif>0))+0.1*max(temp_yield(temp_sif>0));
   plot(temp_yield(temp_sif>0),temp_gpp(temp_sif>0),'ro',...
       start:0.0001:stop,polyval(pp5,start:0.0001:stop),'r-');
   axis([0 0.001 0 50]);
   
   xlabel('SIFyield','FontSize',15,'FontName','Whitney-Book');   
   ylabel('GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   info = [num2str(gpp_doy(doy_i),'%d'),'  ','R2=',num2str(rsq2,' %4.2f'),'  ','RMSE=',num2str(rmse2,' %4.2f')];
   text(pos(1),pos(2),info,'FontSize',20,'FontName','Whitney-Book');
   pos2 = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.15*ht];
   eq = ['GPP=' num2str(pp5(2),'%6.4f') '+' num2str(pp5(1),'%6.4f') '*SIFyield   ' 'p=' num2str(pvalue2,'%6.4f') ];
   text(pos2(1),pos2(2),eq ,'FontSize',15,'FontName','Whitney-Book');
   %=================
   
   fig_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\sif_gpp\',num2str(gpp_doy(doy_i)),'_f760_comp.jpg'));
   fig2_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\sif_gpp_linear\',num2str(gpp_doy(doy_i)),'_f760_linear.jpg'));
   fig3_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\sifyield\',num2str(gpp_doy(doy_i)),'_f760_yield.jpg'));
   fig4_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\sif_gpp_mm\',num2str(gpp_doy(doy_i)),'_f760_mm.jpg'));
   fig5_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\gpp_0328\sifyield_linear\',num2str(gpp_doy(doy_i)),'_f760_yield_linear.jpg'));
   
   saveas(fig,fig_fp);
   close(fig);
   saveas(fig2,fig2_fp);
   close(fig2);
   saveas(fig3,fig3_fp);
   close(fig3);
   saveas(fig4,fig4_fp);
   close(fig4);
   saveas(fig5,fig5_fp);
   close(fig5);
 end
 
 [nelements,centers] = hist(rsq_total(~isnan(rsq_total)),10);
 save('diurnal_figures.mat');
 clearvars
 
 
 
 