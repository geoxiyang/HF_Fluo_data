%% Scratch for making plots

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

load('SIF760daily.mat', 'SIF_mean')
load('hf_barn_2013_env.mat', 'apar_daily')
load('hf_barn_2013_env.mat', 'gpp_day')
load('SIF760daily.mat', 'SIF_1330')
ratio3 = SIF_mean(:,2)./gpp_day(:,1);
yy = ratio3(ratio3(:,1)>0 & ratio3(:,1)<0.15,1);
xx = apar_daily(ratio3(:,1)>0 & ratio3(:,1)<0.15,1);

slm     = slmengine(xx,yy,'knots',3);
model   = slm;

plot(xx,yy,'ko','MarkerSize',12);
[r,m,b] = regression(xx',yy');

hold on

plot([1,max(xx)],[m+b,m*max(xx)+b],'b-');

if strcmpi(model.form,'slm')
  xrange = model.knots([1,end]);
else
  xrange = model.breaks([1,end]);
end

xev = linspace(xrange(1),xrange(2),1001);

% evaluate
if strcmpi(model.form,'slm')
  ypred = slmeval(xev,model);
else
  ypred = ppval(model,xev);
end

% plot the curve
h = plot(xev,ypred);
set(h,'Marker','none', ...
  'Color','r', ...
  'LineStyle','-', ...
  'LineWidth',1);

xlabel('APAR(umol/m^{2}/sec)','FontName','Whitney','FontSize',20);
ylabel('SIF/GPP','FontName','Whitney','FontSize',20);
set(gca,'FontName','Whitney','FontSize',16);

hold off



