%% Read MODIS EOS-HDF file and extract data

import matlab.io.hdfeos.*

% Information to modify

site_lon = -72.1715;
site_lat = 42.5378;

filedir  = '/Volumes/XiYangResearch/Data/2.SatelliteData/2.MODIS/MCD15A2/';
filename = dir('/Volumes/XiYangResearch/Data/2.SatelliteData/2.MODIS/MCD15A2/*.hdf');
gridname = 'MOD_Grid_MOD15A2';

% Define array
dateobs              = nan(length(filename),1);
fieldvalue1          = nan(length(filename),1);
fieldvalue2          = nan(length(filename),1);

% Read files
for i = 1:length(filename)
    
    dateobs(i)       = str2num(filename(i).name(14:16)) * 1.0;
    
    fileid           = matlab.io.hdfeos.gd.open(fullfile(filedir,filename(i).name), 'rdonly');
    gridid           = matlab.io.hdfeos.gd.attach(fileid, gridname);

    dfield1          = 'Fpar_1km';
    dfield2          = 'FparStdDev_1km';

    % Read data from the HDF-EOS2 Grid data field.
    % [data, lon, lat] = matlab.io.hdfeos.gd.readField(gridid, dfield);
    [row , col]      = matlab.io.hdfeos.gd.getPixels(gridid,site_lat,site_lon);
    fieldvalue1(i)   = matlab.io.hdfeos.gd.getPixValues(gridid,row,col,dfield1);
    fieldvalue2(i)   = matlab.io.hdfeos.gd.getPixValues(gridid,row,col,dfield2);

end

final_result = [dateobs,fieldvalue1,fieldvalue2];
final_result(:,1) = final_result(:,1) * 1.0;
final_result(:,2) = final_result(:,2) * 0.01;
final_result(:,3) = final_result(:,3) * 0.01;








% Detach from the grid object.
% matlab.io.hdfeos.gd.detach(gridid);
% 
% Close the HDF-EOS2 grid data file.
% matlab.io.hdfeos.gd.close(fileid);



