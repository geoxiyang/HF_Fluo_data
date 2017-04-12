%% sun_position_site.m
%  Calculate the sun position of a specific site
%  Here we calculate the sun position from June 19 to October 27

clear variables
clc

location.longitude = -72.1899; 
location.latitude = 42.5353; 
location.altitude = 200;

time.year = 2013;
cal_i     = 1;

sza       = nan(1,1);
azi       = nan(1,1);
doy       = nan(1,1);


for day_i = 1:130
   
    tmp_str = strsplit(datestr(datenum(2013,6,20)+day_i-1,'mm/dd/yyyy'),'/');
    time.year = str2double(tmp_str(3));
    time.month = str2double(tmp_str(1));
    time.day = str2double(tmp_str(2));
    
    for hour_i = 1:24
        
        time.hour = -1+hour_i;
        time.min = 30;
        time.sec = 0;
        time.UTC = -5;
        
        sun = sun_position(time, location);
        
        sza(cal_i,1) = sun.zenith;
        azi(cal_i,1) = sun.azimuth;
        doy(cal_i,1) = datenum(str2double(tmp_str(3)),str2double(tmp_str(1)),str2double(tmp_str(2)),hour_i-1,30,0) - datenum(str2double(tmp_str(3)),1,1,0,0,0);
        
        cal_i = cal_i+1;
    end
    
end


