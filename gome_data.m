% process GOME-2 data
% Extract Harvard Forest SIF data (SIF 760)

filepath = 'F:\01.ResearchProject\9.Fluorescence\8.SatelliteData\GOME2\lvl3\*.nc';
filename = dir(filepath);

SIF740_hf = nan(numel(filename),1);
SIF740_std_hf = nan(numel(filename),1);
Par_SIF740_hf = nan(numel(filename),1);
PAR_SIF740_std_hf = nan(numel(filename),1);
cos_SZA_hf = nan(numel(filename),1);

for i = 1:numel(filename) %numel(filename)
    display(['Reading ' filename(i).name]);
    id=netcdf.open(filename(i).name,'NOWRITE');

    varSIF740 = netcdf.inqVarID(id,'SIF_740');
    varSIF740_std = netcdf.inqVarID(id,'SIF_740_std');
    varPar_SIF740 = netcdf.inqVarID(id,'Par_normalized_SIF_740');
    varPar_SIF740_std = netcdf.inqVarID(id,'Par_normalized_SIF_740_std');
    varcos_SZA = netcdf.inqVarID(id,'cos(SZA)');
    varlat = netcdf.inqVarID(id,'latitude');
    varlon = netcdf.inqVarID(id,'longitude');
    
    SIF740 = netcdf.getVar(id,varSIF740);
    SIF740_std = netcdf.getVar(id,varSIF740_std);
    Par_SIF740 = netcdf.getVar(id,varPar_SIF740);
    Par_SIF740_std = netcdf.getVar(id,varPar_SIF740_std);
    cos_SZA = netcdf.getVar(id,varcos_SZA);
    lat = netcdf.getVar(id,varlat);
    lon = netcdf.getVar(id,varlon);
    
    hf_loc = [42.5353,-72.1899];
    lat_sub = knnsearch(lat,hf_loc(1),'K',1);
    lon_sub = knnsearch(lon,hf_loc(2),'K',1);
    
    SIF740_hf(i) = SIF740(lon_sub,lat_sub);
    SIF740_std_hf(i) = SIF740_std(lon_sub,lat_sub);
    Par_SIF740_hf(i) = Par_SIF740(lon_sub,lat_sub);
    PAR_SIF740_std_hf(i) = Par_SIF740_std(lon_sub,lat_sub);
    cos_SZA_hf(i) = cos_SZA(lon_sub,lat_sub); 
    
end

clearvars -except SIF740_hf SIF740_std_hf ...
                  Par_SIF740_hf PAR_SIF740_std_hf ...
                  cos_SZA_hf