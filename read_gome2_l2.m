function [lat, lon, sza, vza, refl670, refl780, resid_rms, f_c, sif_740, wav] = read_gome2_l2(filename)
%% Code processess chlorophyll fluorescence data (netcdf) downloaded from NASA GOME-2 data site provided by J. Joiner
% website: http://acdb-ext.gsfc.nasa.gov/People/Joiner/my_gifs/GOME_F/GOME-F.htm
% Adapted from IDL code provided from website listed above
% Y.P. Shiga 7/29/13 yshiga@stanford.edu
%%
display(['Reading ' filename]);
id=netcdf.open(filename,'NOWRITE');   % netcdf open file command
varwav=netcdf.inqVarID(id,'Wavelengths');   % sets "Wavelengths" variable id name
if varwav>=0
    % sets variable id names
    varFsc = netcdf.inqVarID(id,'SIF_740'); 
    varFsc_std=netcdf.inqVarID(id,'cloud_fraction');
    varrf_avg=netcdf.inqVarID(id,'refl670');
    varrf_std=netcdf.inqVarID(id,'refl780');
    varcossza=netcdf.inqVarID(id,'SZA');
    var=netcdf.inqVarID(id,'VZA');
    varresid=netcdf.inqVarID(id,'residual');
    varlat = netcdf.inqVarID(id,'latitude');
    varlon = netcdf.inqVarID(id,'longitude');
    %reads data (from above variable id)
    wav= netcdf.getVar(id,varwav);
    sif_740 = netcdf.getVar(id,varFsc);
    f_c = netcdf.getVar(id,varFsc_std);
    refl670 = netcdf.getVar(id,varrf_avg);
    refl780 = netcdf.getVar(id,varrf_std);
    sza = netcdf.getVar(id,varcossza);
    vza = netcdf.getVar(id,var);
    resid_rms = netcdf.getVar(id,varresid);
    lat = netcdf.getVar(id,varlat);
    lon = netcdf.getVar(id,varlon);
else display(['Error Reading file/No data ' filename])
end
netcdf.close(id)
end

