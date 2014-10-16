%% cossza corrected SIF


load('SIF760daily.mat','SIF_0930','SIF_mean');
load('hf_barn_2013_env.mat','cossza_day','par_daily','apar_daily');

SIF_mean_sza(:,1)   = SIF_mean(:,1);
for ii = 2:4
    SIF_mean_sza(:,ii) = SIF_mean(:,ii)./cossza_day';
end
save('SIF760daily.mat','SIF_mean_sza','-append');
% DOY = SIF_mean(SIF_mean(:,2)>0,1);
% SIF_mean_sza = SIF_mean_sza(SIF_mean_sza>0);
% 
% [AX,H1,H2] = plotyy(DOY,SIF_mean_sza,170:1:299,apar_daily);


% corr(SIF_mean_sza(SIF_mean(:,2)>0 & ~isnan(apar_daily)),apar_daily(SIF_mean(:,2)>0 & ~isnan(apar_daily)))

