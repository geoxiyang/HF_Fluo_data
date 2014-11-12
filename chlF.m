%calculate ChlF
clear all
clc

tic

coeff = importdata('/Volumes/XiYangResearch/Projects/9.Fluorescence/1.coefficients/coeff.txt');

datapath = '/Volumes/XiYangResearch/Data/13.Fluorescence/2014/measurements/';

rad_path = dir(fullfile(strcat(datapath,'rad/*.txt')));
irrad_path = dir(fullfile(strcat(datapath,'irrad/*.txt')));
file_no = size(rad_path);
time = NaN(file_no);
f760 = NaN(file_no);
rsq = NaN(file_no);
rmse = NaN(file_no);
ref = NaN(file_no(1),2048); % 2048 sampling intervals

for i = 1:file_no(1)  
    
    fname = fullfile(strcat(datapath,'rad/',rad_path(i).name));
    fid = fopen(fname);
    headline = fgetl(fid);
    fclose(fid);
    temp_rad = dlmread(fname, ' ',3, 0);
    
    fname0 = fullfile(strcat(datapath,'irrad/',irrad_path(i).name));
    fid0 = fopen(fname0);
    headline0 = fgetl(fid0);
    fclose(fid0);
    temp_irrad = dlmread(fname0, ' ',3, 0);
    
    raw_time = strread(rad_path(i).name,'%s','delimiter','_');
    raw_time_irrad = strread(irrad_path(i).name,'%s','delimiter','_');
   
    rad_time_tmp    = datenum([str2double(raw_time{2}) str2double(raw_time{3})...
        str2double(raw_time{4}) str2double(raw_time{5}) str2double(raw_time{6}) 0]) - datenum( [2014 1 1 0 0 0]) + 1;
    irrad_time_tmp  = datenum([str2double(raw_time_irrad{2}) str2double(raw_time_irrad{3})...
        str2double(raw_time_irrad{4}) str2double(raw_time_irrad{5}) str2double(raw_time_irrad{6}) 0]) - datenum( [2014 1 1 0 0 0]) + 1;
    
    tenmin = datenum( [2014 1 1 0 10 0]) - datenum( [2014 1 1 0 0 0]);
    if (rad_time_tmp - irrad_time_tmp > 0) || (rad_time_tmp - irrad_time_tmp > tenmin)
        continue
    end
    time(i) = datenum([str2double(raw_time{2}) str2double(raw_time{3})...
         str2double(raw_time{4}) str2double(raw_time{5}) str2double(raw_time{6}) 0]) - datenum( [2014 1 1 0 0 0])+1;
    
    %Reflectance
    ref(i,:) = (temp_rad(:,3).*coeff(:,2))./(temp_irrad(:,3).*coeff(:,1));
    wl       = temp_rad(:,5);
    
    
    
    %For O2A: 759-767.76, centered 760;
    roi = temp_rad(:,5)>759.00 & temp_rad(:,5)<767.76;
   
    wavelength = temp_rad(roi,5);
    rad = temp_rad(roi,3).*coeff(roi,2);
    irrad = temp_irrad(roi,3).*coeff(roi,1);
    xdata = [wavelength irrad]';
    ydata = rad';
    
    % SFM fitting 
    x0 = [0.1,0.1,0.1,0.1];
    lb = [-inf,-inf,-inf,-inf];
    ub = [inf,inf,inf,inf];
    options = optimset('TolX',1e-10,'TolFun',1e-10,'Display','off');
    [x,SSresid] = lsqcurvefit(@radtrans,x0,xdata,ydata,lb,ub,options);
    
    SStotal = (length(ydata)-1) * var(ydata);
    rsq(i) = 1 - SSresid/SStotal;
    rmse(i) = (SSresid/(length(ydata)))^(0.5);
    %For O2A: center 760
    f760(i) = 760.0*x(4)+x(3);
    
end

% Store reflectance results
ref_final       = [time, ref];
ref_final_time  = sortrows(ref_final,1);
raw_ref_result  = ref_final_time;

% save('SIF760_result.mat','raw_ref_result','wl','-append');

% Store SIF results
final_result = [time, f760, rsq, rmse];

% Time selection
final_result_time = sortrows(final_result,1);

raw_final_result = final_result_time;
 