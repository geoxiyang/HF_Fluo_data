%% make figures for HF fluorescence paper

% Read GOME-2 site data

load('SIF760daily.mat','SIF_0930');
load('hf_barn_2013_env.mat','cossza_0930');

% SZA correction
% SIF_0930(:,2) = SIF_0930(:,2)./cossza_0930';
% SIF_0930(:,3) = SIF_0930(:,3)./cossza_0930';
% SIF_0930(:,4) = SIF_0930(:,4)./cossza_0930';

SIF_0930 = SIF_0930(~ismember(SIF_0930(:,1),[172,205,187,189,261]),:);

SIF_0930_1 = SIF_0930(SIF_0930(:,2)>0.0,1:2);
SIF_0930_2 = SIF_0930(SIF_0930(:,3)>0.0,[1,3]);

ofigure  = '/Users/xiyang/Dropbox/Mypaper/6.2014-Fluorescence/JPG/';
filename = '/Volumes/XiYangResearch/Projects/14.CropSIF/result/weekly_v25/HF_3x3.dat';
rawdata = importdata(filename);

zero_pad = zeros(2,10);
rawdata = vertcat(zero_pad,rawdata);

year = rawdata(:,1);
month = rawdata(:,2);
doy = rawdata(:,4);
sif = rawdata(:,5);
sif_sd = rawdata(:,6);
sif_par = rawdata(:,7);
sif_par_sd = rawdata(:,8);
ndvi = rawdata(:,9);
ndvi_sd = rawdata(:,10);

[row col] = size(rawdata);

% Change the reshape size if processing other sizes of data
sif_year = reshape(sif,row/7,7);
ndvi_year = reshape(ndvi,row/7,7);
sif_par_year = reshape(sif_par,row/7,7);

sif_sd_year = reshape(sif_sd,row/7,7);
ndvi_sd_year = reshape(ndvi_sd,row/7,7);
sif_par_sd_year = reshape(sif_par_sd,row/7,7);

%calculate the mean and sd with 2012
sif_all_sg = zeros(row/7,7); 
for ii = 1: 7
   sif_all_sg(:,ii) =  sgolayfilt(sif_year(:,ii),3,7);   %change to sif_par
end
mean_sif_all_sg = mean(sif_all_sg,2);
sd_sif_all = std(sif_all_sg,0,2);


for jj = ceil(170/7):(row/7)
    SIF_0930_SUN_TEMP = SIF_0930_2(SIF_0930_2(:,1) >= ((jj-1)*7+1) &...
                                   SIF_0930_2(:,1) <  (jj*7+1),2);
    SIF_0930_SUN_WEEK(jj) = median(SIF_0930_SUN_TEMP);
end
set(gca,'xlim',[170 299]);

% [h1,hp] = boundedline(1:7:row, mean_sif_all_sg.*0.582, sd_sif_all, 'ro-','transparency', 0.1);

[h1,hp] = boundedline(1:7:row, sif_year(:,7).*0.582, sd_sif_all, 'ro-','transparency', 0.1);
set(h1,'MarkerSize',16,'LineWidth',2);

hold on
h3 = plot(1:7:row,SIF_0930_SUN_WEEK,'bo','MarkerSize',16,'LineWidth',2);

hold off
%plot(SIF_0930_SUN_WEEK(ceil(170/7):(row/7)),mean_sif_all_sg(ceil(170/7):(row/7)),'b*')

ylabel('SIF(mw/m^{2}/sr/nm)','FontSize',16,'FontName','Whitney-Book');
xlabel('Day of Year','FontSize',16,'FontName','Whitney-Book');
set(gca,'FontSize',16,'FontName','Whitney-Book');
legend([h1,h3],{'GOME-2 SIF','Ground SIF'});

t_1 = SIF_0930_SUN_WEEK(ceil(170/7):(row/7))';
t_2 = sif_year(ceil(170/7):(row/7),7);

corr(t_1(~isnan(t_1)),t_2(~isnan(t_1)))


print('-dpng','-r300',[ofigure 'GOME-2-Ground.png']);


% 
% t_1 = SIF_0930_SUN_WEEK(ceil(170/7):(row/7))';
% t_2 = sif_all_sg(ceil(170/7):(row/7),7);
% 
% corr(t_1(~isnan(t_1)),t_2(~isnan(t_1)))
% 
% plot(t_1(~isnan(t_1)),t_2(~isnan(t_1)),'ko')
% ylabel('GOME-2 SIF(mw/m2/sr/nm)','FontSize',16,'FontName','Whitney-Book');
% xlabel('Ground SIF(mw/m2/sr/nm)','FontSize',16,'FontName','Whitney-Book');
% set(gca,'FontSize',16,'FontName','Whitney-Book');




