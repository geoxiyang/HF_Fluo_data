%% Load fluorescence data and plot them

%% 1. Load data

work_dir    = '/Volumes/XiYangResearch/src/SCOPE_v1.40beta/SCOPE_v1.40/output/HF_2013_2014-08-14-0940/';
wl          = dlmread([work_dir 'wl.dat'],'',2,0);
fEnergy     = dlmread([work_dir 'fluorescence.dat'],'',2,0);

load('/Volumes/XiYangResearch/src/HF_Fluo_data/SIF760daily.mat','halfhourly_result');

ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/figure.png';


%% 2. Plot data

wl          = wl(1,wl(1,:) >=640 & wl(1,:)<=850);
scope_f760  = fEnergy(:,wl(1,:) == 760);

day         = 175;
measure_f760= halfhourly_result(floor(halfhourly_result(:,1)) == day,:);

scope_f760 = [measure_f760(:,1),scope_f760];
plot(scope_f760(:,1)-floor(scope_f760(:,1)),scope_f760(:,2),'bo',measure_f760(:,1)-floor(measure_f760(:,1)),measure_f760(:,2),'r*')

xlabel('Hour','FontSize',20,'FontName','Helvetica');
ylabel('SIF(w m-2 nm-1 sr-1)','FontSize',20,'FontName','Helvetica');
set(gca,'FontSize',20,'FontName','Helvetica')
legend('Modelled SIF','Measured SIF');

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', ofigure);
