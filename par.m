% for i = 1:128
%     doy_fs_sub = find(doy2 == doy1(i)); 
%    par_fs(i,1) = par(doy_fs_sub);  
%     
% end

% par15fill = zeros(15681,1);
% 
% for i=1:15681
%    if isnan(final_result_time(i,1))
%         par15fill(i,1) = NaN;         
%         continue
%    end
%     
%    [subs,dists] = knnsearch(par15doy,final_result_time(i,1),'K',2);
%    
%    weight1 = 1/abs(par15doy(subs(1))-final_result_time(i,1));
%    weight2 = 1/abs(par15doy(subs(2))-final_result_time(i,1));
%    
%    par15fill(i,1) = (weight1*par15min(subs(1))+weight2*par15min(subs(2)))/(weight1+weight2);
%     
%     
%     
% end
% save('SIF760_result.mat');

load('SIF760_result.mat');
load('SIF760_envi.mat','par15fill','vpd15fill','airt15fill');

 %doy = [184,199,215,216,217,218,221,222,226,227,231,232,233,246,248,249,253,262,265,267,273,275];
%doy = [186,192,193,198,200,201,204,213,224,268,271,272,274,289,293];
doy = double(uni_doy(:,1)');
nonzero_ind = find(doy);
doy = doy(nonzero_ind);
doyn = size(doy);

% testa = [time, f760, rsq,rmse];
% raw_final_result = sortrows(testa,1);

rsq_total = nan(numel(doy),1);
p_total = nan(numel(doy),1);
rmse_total = nan(numel(doy),1);

for doy_i=4:4 %1:doyn(2)


   lb = doy(doy_i)+0.;
   ub = doy(doy_i)+1.;
   temp_time = 24.0*(raw_final_result(raw_final_result(:,1) >= lb & raw_final_result(:,1) <= ub & raw_final_result(:,3) >=0.90, 1) - doy(doy_i));
   temp_array = raw_final_result(raw_final_result(:,1) >= lb & raw_final_result(:,1) <= ub & raw_final_result(:,3) >=0.90, 2);
   temp_par15 = par15fill(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub & raw_final_result(:,3) >=0.90, 1);
   temp_yield = temp_array/temp_par15;
   
   
   size_test = size(temp_time);
   if size_test < 5
       continue
   end
%    temp_vpd15 = vpd15fill(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 1);
%    temp_airt15 = airt15fill(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 1);
%    temp_gpp = gppfill(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 1);
%    
%    plot(temp_time,temp_array,'-o');
   
   
   %temp_time,temp_par15,temp_time,temp_vpd15,temp_time,temp_airt15);
   

   %Figure 1   
   fig = figure;
   set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 10 6]);
   set(fig,'Visible','Off');
   [AX,H1,H2] = plotyy(temp_time,temp_array,temp_time,temp_par15);
%    [AX,H1,H2] = plotyy(temp_time,temp_array,temp_time,temp_gpp);
   set(H1,'marker','O');
   set(H2,'marker','O');
   set(get(AX(1),'Ylabel'),'String','SIF(mw/m2/sr/nm)','FontSize',15,'FontName','Whitney-Book');
   set(get(AX(2),'Ylabel'),'String','PAR(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
%    set(get(AX(2),'Ylabel'),'String','GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
   set(AX(1),'YLim',[0 5],'FontSize',12,'FontName','Whitney-Book');
   set(AX(2),'YLim',[0 2000],'FontSize',12,'FontName','Whitney-Book');
%    set(AX(2),'YLim',[0 50],'FontSize',12,'FontName','Whitney-Book');

   xlabel('Hours','FontSize',15,'FontName','Whitney-Book');
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   text(pos(1),pos(2),num2str(doy(doy_i)),'FontSize',20,'FontName','Whitney-Book');
   %=================
   
   %Figure 2
   sub1 = ~isnan(temp_par15);
%    sub1 = ~isnan(temp_gpp);
   sub2 = ~isnan(temp_array);
   sub = sub1 & sub2;
   
   fitresult = linearfit(temp_par15(sub),temp_array(sub));
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
   set(fig2,'Visible','Off');
   plot(temp_par15,temp_array,'ro',0:1:2000,polyval([slope,intercept],0:1:2000),'r-');
%    plot(temp_array,temp_gpp,'ro',0:1:50,polyval([slope,intercept],0:1:50),'r-');
   axis([0 2000 0 5]);
%    axis([0 50 0 5]);
   
   xlabel('PAR(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');
%    xlabel('GPP(umol/m2/second)','FontSize',15,'FontName','Whitney-Book');

   ylabel('SIF(mw/sr/m2/nm)','FontSize',15,'FontName','Whitney-Book');   
   
   xlim_fig = xlim;
   ylim_fig = ylim;
   
   wdth = xlim_fig(2)-xlim_fig(1);
   ht = ylim_fig(2)-ylim_fig(1);
   pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.05*ht];
   info = [num2str(doy(doy_i),'%d'),'  ','R2=',num2str(rsq,' %4.2f')];
   text(pos(1),pos(2),info,'FontSize',20,'FontName','Whitney-Book');
   pos2 = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.15*ht];
   eq = ['y=' num2str(intercept,'%6.4f') '+' num2str(slope,'%6.4f') '*x   ' 'p=' num2str(pvalue,'%6.4f') ];
   text(pos2(1),pos2(2),eq ,'FontSize',15,'FontName','Whitney-Book');
   %=================

   
   fig_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\rsq_control\',num2str(doy(doy_i)),'_f760_comp.jpg'));
   fig2_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\rsq_control\scatter\',num2str(doy(doy_i)),'_f760_scatter.jpg'));
   
   saveas(fig,fig_fp);
   close(fig);
   saveas(fig2,fig2_fp);
   close(fig2);
end
 

save('SIF760_result.mat','-append');
% for uni_i = 1:doy_size(1)  %1:doy_size(1)
%     
%    lb = double(uni_doy(uni_i))+0.;
%    ub = double(uni_doy(uni_i))+1.;
%    temp_time = final_result_time(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 1);
%    temp_array = final_result_time(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 2);
%    temp_par15 = par15fill(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 1);
%    
%    
%    fig = figure;
%    set(fig,'Visible','off');
%    ylim([0 5]);
%    plot(temp_par15,temp_array,'r.','MarkerSize',20);
%    
%    
% %    fig_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\',num2str(uni_doy(uni_i)),'_f687.jpg'));
%     fig_fp = fullfile(strcat('F:\01.ResearchProject\9.Fluorescence\4.JPG\',num2str(uni_doy(uni_i)),'_PAR_SIF.jpg'));
% % 
%     saveas(fig,fig_fp);
%     close(fig);
%     
% %    lb = double(uni_doy(uni_i))+0.;
% %    ub = double(uni_doy(uni_i))+1.;
% %    
% %    mid_day_result(1,uni_i) = double(uni_doy(uni_i));
% %    temp_array = final_result_time(final_result_time(:,1) >= lb & final_result_time(:,1) <= ub, 2);
% %    if isempty(temp_array) 
% %         mid_day_result(2,uni_i) = NaN;
% %    else
% %         mid_day_result(2,uni_i) = nanmax(temp_array);
% %    end
% end