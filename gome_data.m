% process GOME-2 data
% Extract Harvard Forest SIF data (SIF 760)

clear variables
clc

filepath = '/Volumes/XiYangResearch/Data/2.SatelliteData/1.GOME-2/GOME_F/GOME_F/MetOp-A/level3/';%'/Volumes/XiYangResearch/Data/2.SatelliteData/1.GOME-2/lvl3_v26_monthly/';
cd(filepath)
filename = dir([filepath '*.nc']);

% SIF740_hf = nan(numel(filename),);
% SIF740_std_hf = nan(numel(filename),1);
% Par_SIF740_hf = nan(numel(filename),1);
% PAR_SIF740_std_hf = nan(numel(filename),1);
% cos_SZA_hf = nan(numel(filename),1);

SIF740          = nan(numel(filename),720,360);
SIF740_SD       = nan(numel(filename),720,360);
Par_SIF740      = nan(numel(filename),720,360);
Par_SIF740_SD   = nan(numel(filename),720,360);
cos_SZA         = nan(numel(filename),720,360);
timeym          = nan(numel(filename),2);
counts_gome2    = nan(numel(filename),720,360);

for i = 1:numel(filename) %numel(filename)
    display(['Reading ' filename(i).name]);
    id=netcdf.open(filename(i).name,'NOWRITE');

    varSIF740           = netcdf.inqVarID(id,'SIF_740');
    varSIF740_std       = netcdf.inqVarID(id,'SIF_740_std');
    varPar_SIF740       = netcdf.inqVarID(id,'Par_normalized_SIF_740');
    varPar_SIF740_std   = netcdf.inqVarID(id,'Par_normalized_SIF_740_std');
    varcos_SZA          = netcdf.inqVarID(id,'cos(SZA)');
    varlat              = netcdf.inqVarID(id,'latitude');
    varlon              = netcdf.inqVarID(id,'longitude');
    varCounts           = netcdf.inqVarID(id,'Counts');
    
    SIF740_raw          = netcdf.getVar(id,varSIF740);
    SIF740_std_raw      = netcdf.getVar(id,varSIF740_std);
    Par_SIF740_raw      = netcdf.getVar(id,varPar_SIF740);
    Par_SIF740_std_raw  = netcdf.getVar(id,varPar_SIF740_std);
    cos_SZA_raw         = netcdf.getVar(id,varcos_SZA);
    counts_gome2_raw    = netcdf.getVar(id,varCounts);
    
    
    lat                 = netcdf.getVar(id,varlat);
    lon                 = netcdf.getVar(id,varlon);
    year                = str2double(filename(i).name(end-13:end-10));
    month               = str2double(filename(i).name(end-09:end-08));

% This is for a single point extraction
%     hf_loc = [42.5353,-72.1899];
%     lat_sub = knnsearch(lat,hf_loc(1),'K',1);
%     lon_sub = knnsearch(lon,hf_loc(2),'K',1);
%     
%     SIF740_hf(i) = SIF740_raw(lon_sub,lat_sub);
%     SIF740_std_hf(i) = SIF740_std_raw(lon_sub,lat_sub);
%     Par_SIF740_hf(i) = Par_SIF740_raw(lon_sub,lat_sub);
%     PAR_SIF740_std_hf(i) = Par_SIF740_std_raw(lon_sub,lat_sub);
%     cos_SZA_hf(i) = cos_SZA_raw(lon_sub,lat_sub); 
    

% Pack all monthly SIF in a single file
% Also store year/month data in a variable

    timeym(i,1)                 = year;
    timeym(i,2)                 = month;
    SIF740(i,:,:)               = SIF740_raw;
    SIF740_SD(i,:,:)            = SIF740_std_raw;
    Par_SIF740(i,:,:)           = Par_SIF740_raw;
    Par_SIF740_SD(i,:,:)        = Par_SIF740_std_raw;
    cos_SZA(i,:,:)              = cos_SZA_raw;
    counts_gome2(i,:,:)         = counts_gome2_raw;
    
end

save('/Volumes/XiYangResearch/Projects/4.DiurnalLUE/2.Matlab/gome_monthly_v26_MetOpA.mat')