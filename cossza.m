%% calculate cos(SZA) - consine of sun zenith angle

year        = 2013;
day         = 170:1:299;
hours       = 10:1:16;

latitude    = 42.54;
longitude   = -72.17;
altitude    = 200;

% Input data
location.longitude      = longitude; 
location.latitude       = latitude; 
location.altitude       = altitude;
time.year               = year;
time.min                = 30;
time.sec                = 0;
time.UTC                = -5;

cossza_day  = zeros(1,numel(day));
cossza_0930 = zeros(1,numel(day));

for ii = 1:numel(day)

[YY,MM,DD,HH] = datevec(datenum(year,1,day(ii)));    
time.month              = MM;
time.day                = DD;

time.hour = 9;
sun = sun_position(time, location);
cossza_0930(1,ii) = cosd(sun.zenith);

% cossza_hour = zeros(1,numel(hours));
%     for jj = 1:numel(hours)
%         time.hour = hours(jj);
%         sun = sun_position(time, location);
%         cossza_hour(1,jj) = cosd(sun.zenith);
%     end
%cossza_day(ii) = mean(cossza_hour);

end

plot(cossza_0930,'ko');




save('hf_barn_2013_env.mat','cossza_0930','-append')



