%% Scratch for making plots

load('HF_2013_GPP.mat', 'gpp_day')
load('SIF760daily.mat', 'halfhourly_result','SIF_mean')
SIF_cube = reshape(halfhourly_result(:,2),48,130);
SIF_cube = SIF_cube';
SIF_cube(SIF_cube<0) = NaN;
SIF_cube(SIF_cube<=0) = NaN;
sif_gome_local = nanmean(SIF_cube(:,18:20),2);
sif_gosat_local = nanmean(SIF_cube(:,26:28),2);


sub5 = sif_gosat_local > 0 & gpp_day > 0 & sif_gosat_local < 2.5;
sub4 = sif_gome_local > 0 & gpp_day > 0;
sub6 = SIF_mean(:,2) > 0 & gpp_day >0 & SIF_mean(:,2) < 2.5;

plot(SIF_mean(sub6,2),gpp_day(sub6),'k.','MarkerSize',24)
hold on
plot(sif_gosat_local(sub5),gpp_day(sub5),'m.','MarkerSize',24)
plot(sif_gome_local(sub4),gpp_day(sub4),'g.','MarkerSize',24)

ylim([0 18]);

% SIF_mean
xint = (nanmax(SIF_mean(sub6,2)) - nanmin(SIF_mean(sub6,2)))/100.0;
xfit = nanmin(SIF_mean(sub6,2)):xint:nanmax(SIF_mean(sub6,2));
x    = [ones(numel(SIF_mean(sub6,2)),1),SIF_mean(sub6,2)];
[b,bint,r,rint,stats] = regress(gpp_day(sub6),x);
yfit = b(1) + b(2) * xfit;
plot(xfit,yfit,'k-','LineWidth',2);
pos1 = [0.1,17.2];
text(pos1(1),pos1(2),['GPP=',num2str(b(1),'% 5.2f'),'+',num2str(b(2),'% 5.2f'),'\times','SIF',...
                       '  ',...
                       'r^{2}=',num2str(corr(SIF_mean(sub6,2),gpp_day(sub6))^2,'% 5.2f'),...
                       '  ',...
                       'p=',num2str((stats(3)),'% 6.4f')],'FontSize',20,'FontName','Whitney','Color','k');

% SIF_gosat_local
xint = (nanmax(sif_gosat_local(sub5)) - nanmin(sif_gosat_local(sub5)))/100.0;
xfit = nanmin(sif_gosat_local(sub5)):xint:nanmax(sif_gosat_local(sub5));
x    = [ones(numel(sif_gosat_local(sub5)),1),sif_gosat_local(sub5)];
[b,bint,r,rint,stats] = regress(gpp_day(sub5),x);
yfit = b(1) + b(2) * xfit;
plot(xfit,yfit,'m-','LineWidth',2);
pos2 = [0.1,15.9];
text(pos2(1),pos2(2),['GPP=',num2str(b(1),'% 5.2f'),'+',num2str(b(2),'% 5.2f'),'\times','SIF',...
                      '  ',...
                      'r^{2}=',num2str(corr(sif_gosat_local(sub5),gpp_day(sub5))^2,'% 5.2f'),...
                      '  ',...
                      'p=',num2str((stats(3)),'% 6.4f')],'FontSize',20,'FontName','Whitney','Color','m');

% SIF_gome_local

xint = (nanmax(sif_gome_local(sub4)) - nanmin(sif_gome_local(sub4)))/100.0;
xfit = nanmin(sif_gome_local(sub4)):xint:nanmax(sif_gome_local(sub4));
x    = [ones(numel(sif_gome_local(sub4)),1),sif_gome_local(sub4)];
[b,bint,r,rint,stats] = regress(gpp_day(sub4),x);
yfit = b(1) + b(2) * xfit;
plot(xfit,yfit,'g-','LineWidth',2);
pos3 = [0.1,14.6];
text(pos3(1),pos3(2),['GPP=',num2str(b(1),'% 5.2f'),'+',num2str(b(2),'% 5.2f'),'\times','SIF',...
                      '  ',...
                      'r^{2}=',num2str(corr(sif_gome_local(sub4),gpp_day(sub4))^2,'% 5.2f'),...
                      '  ',...
                      'p=',num2str((stats(3)),'% 6.4f')],'FontSize',20,'FontName','Whitney','Color','g');


xlabel('SIF(mw/m^{2}/sr/nm)','FontName','Whitney','FontSize',24);
ylabel('GPP(g C/m^{2}/day)','FontName','Whitney','FontSize',24);
set(gca,'FontName','Whitney','FontSize',20);
title('');

hold off

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', '/Users/xiyang/Dropbox/Mypaper/6.2014-Fluorescence/resub2/figures/SIF_GPP_multiple_times.png');
% close(gcf);







% SIF_0930 vs. GPP_0930
% plot(SIF_0930(SIF_0930(:,3)>0,3),gpp_day(SIF_0930(:,3)>0,1),'ko','MarkerSize',12);
% 
% xlabel('SIF@09:30','FontName','Whitney','FontSize',20);
% ylabel('GPP@09:30','FontName','Whitney','FontSize',20);
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIF0930_GPP0930.png');
% close(gcf);
% 
% % SIF_1400 vs. GPP_1400
% plot(SIF_1400(SIF_1400(:,3)>0,3),gpp_day(SIF_1400(:,3)>0,1),'ko','MarkerSize',12);
% 
% xlabel('SIF@14:00','FontName','Whitney','FontSize',20);
% ylabel('GPP@14:00','FontName','Whitney','FontSize',20);
% set(gca,'FontName','Whitney','FontSize',16);
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIF1400_GPP1400.png');
% close(gcf);

% load('SIF760daily.mat', 'SIF_mean')
% load('hf_barn_2013_env.mat', 'apar_daily')
% load('hf_barn_2013_env.mat', 'gpp_day')
% load('SIF760daily.mat', 'SIF_1330')
% ratio3 = SIF_mean(:,2)./gpp_day(:,1);
% yy = ratio3(ratio3(:,1)>0 & ratio3(:,1)<0.15,1);
% xx = apar_daily(ratio3(:,1)>0 & ratio3(:,1)<0.15,1);
% 
% slm     = slmengine(xx,yy,'knots',3);
% model   = slm;
% 
% plot(xx,yy,'ko','MarkerSize',12);
% [r,m,b] = regression(xx',yy');
% 
% hold on
% 
% plot([1,max(xx)],[m+b,m*max(xx)+b],'b-');
% 
% if strcmpi(model.form,'slm')
%   xrange = model.knots([1,end]);
% else
%   xrange = model.breaks([1,end]);
% end
% 
% xev = linspace(xrange(1),xrange(2),1001);
% 
% % evaluate
% if strcmpi(model.form,'slm')
%   ypred = slmeval(xev,model);
% else
%   ypred = ppval(model,xev);
% end
% 
% % plot the curve
% h = plot(xev,ypred);
% set(h,'Marker','none', ...
%   'Color','r', ...
%   'LineStyle','-', ...
%   'LineWidth',1);
% 
% xlabel('APAR(umol/m^{2}/sec)','FontName','Whitney','FontSize',20);
% ylabel('SIF/GPP','FontName','Whitney','FontSize',20);
% set(gca,'FontName','Whitney','FontSize',16);
% 
% hold off



