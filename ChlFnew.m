%% Calculate SIF at 760nm
%  This is the new version that automatically finds the radiance
%  measurement 10 minutes after the irradiance measurements

%  Calculates SIF using both SFM (Meroni et al., 2009)
%  and SVD (Guanter et al., 2012)


% clear all
%clc

tic

coeff = importdata('/Volumes/XiYangResearch/Projects/9.Fluorescence/1.coefficients/coeff.txt');

datapath = '/Volumes/XiYangResearch/Data/13.Fluorescence/2013/';

rad_path    = dir(fullfile(strcat(datapath,'rad/*.txt')));
irrad_path  = dir(fullfile(strcat(datapath,'irrad/*.txt')));
file_no     = size(rad_path);
time        = NaN(file_no);
f760        = NaN(file_no);
rsq         = NaN(file_no);
rmse        = NaN(file_no);
ref         = NaN(file_no(1),2048); % 2048 sampling intervals
sif_svd     = NaN(file_no);
sif_error   = NaN(file_no);


irrad_time  = NaN(file_no);
rad_time    = NaN(file_no);

% Store the time of each measurement
% for i = 1:file_no(1)
% 	fname   = fullfile(strcat(datapath,'rad/',rad_path(i).name));
%     fname0  = fullfile(strcat(datapath,'irrad/',irrad_path(i).name));
%     
%     raw_time = strread(rad_path(i).name,'%s','delimiter','_');
%     raw_time_irrad = strread(irrad_path(i).name,'%s','delimiter','_');
%    
%     rad_time(i)    = datenum([str2double(raw_time{2}) str2double(raw_time{3})...
%         str2double(raw_time{4}) str2double(raw_time{5}) str2double(raw_time{6}) 0]) - datenum( [2013 1 1 0 0 0]) + 1;
%     irrad_time(i)  = datenum([str2double(raw_time_irrad{2}) str2double(raw_time_irrad{3})...
%         str2double(raw_time_irrad{4}) str2double(raw_time_irrad{5}) str2double(raw_time_irrad{6}) 0]) - datenum( [2013 1 1 0 0 0]) + 1;
% end
% 
% save('hf_2013_svd.mat','rad_time','irrad_time');


load('hf_2013_svd.mat','rad_time','irrad_time')

jj=0;

%4461 = DOY 175 around 13:30

for j = 4461:4461 %1:file_no(1)
	tenmin          = datenum([2013 1 1 0 10 0]) - datenum([2013 1 1 0 0 0]);
    irrad_time_tmp  = irrad_time(j);
    time(j)         = irrad_time(j);
    rad_time_sub    = find((rad_time - irrad_time_tmp)>0 & (rad_time - irrad_time_tmp)<tenmin);

    % if the number of subscript that fit the above criteria is less than
    % one, then quit the loop
    if ~(length(rad_time_sub) == 1)
        continue
    end
    
    % open the irradiance file
    fname_irrad = fullfile(strcat(datapath,'irrad/',irrad_path(j).name));
    temp_irrad  = dlmread(fname_irrad, ' ',3, 0);
    % open the radiance file
    fname_rad   = fullfile(strcat(datapath,'rad/',rad_path(rad_time_sub(1)).name));
    temp_rad    = dlmread(fname_rad, ' ',3, 0);
     
    ref(j,:) = (temp_rad(:,3).*coeff(:,2))./(temp_irrad(:,3).*coeff(:,1));
    wl       = temp_rad(:,5);
    
    %For O2A: 759-767.76, centered 760;
    roi         = temp_rad(:,5)>759.00 & temp_rad(:,5)<767.76;     
    wavelength  = temp_rad(roi,5);
    rad         = temp_rad(roi,3).*coeff(roi,2);
    irrad       = temp_irrad(roi,3).*coeff(roi,1);
    xdata       = [wavelength irrad]';
    ydata       = rad';

    % SFM fitting 
    x0 = [0.1,0.1,0.1,0.1];
    lb = [-inf,-inf,-inf,-inf];
    ub = [inf,inf,inf,inf];
    options = optimset('TolX',1e-10,'TolFun',1e-10,'Display','off');
    [x,SSresid] = lsqcurvefit(@radtrans,x0,xdata,ydata,lb,ub,options);
    
    SStotal = (length(ydata)-1) * var(ydata);
    rsq(j) = 1 - SSresid/SStotal;
    rmse(j) = (SSresid/(length(ydata)))^(0.5);
    %For O2A: center 760
    f760(j) = 760.0*x(4)+x(3);   
    
    spectra = struct('irrad', temp_irrad(:,3)',...
                     'rad', temp_rad(:,3)',...
                     'ircoeff', coeff(:,1),...
                     'rcoeff', coeff(:,2),...
                     'wl', temp_rad(:,5));
    WL_range = [745.00,775.00];
    SIF = SIF_SVD(spectra,WL_range);
    sif_svd(j) = SIF.SIF;
    sif_error(j) = SIF.SIF_error;
%    sif_day = [time_day,sif_svd(rad_idx)];
    
    
    
    
end




% final_svd = [time,sif_svd,sif_error];
% final_svd = sortrows(final_svd,1);
% 
% save('hf_2013_svd.mat','rad_time','irrad_time','final_svd');

% SVD method
% for kk = 175:175
%     
%     irrad_idx = find(irrad_time >= kk+0.3555 & irrad_time <= kk+0.3565);
%     temp_irrad_day = NaN(length(irrad_idx),2048);
%     
%     if isempty(irrad_idx)
%         continue
%     else
%        for irrad_ii = 1:length(irrad_idx) 
%         fname_irrad = fullfile(strcat(datapath,'irrad/',irrad_path(irrad_idx(irrad_ii)).name));
%         temp_irrad  = dlmread(fname_irrad, ' ',3, 0); 
%         temp_irrad_day(irrad_ii,:) = temp_irrad(:,3);
%        end
%     end
%     
%     rad_idx   = find(rad_time >= kk+0.3555 & rad_time <= kk+0.3565);
%     temp_rad_day   = NaN(length(rad_idx),2048);
% 
%     if isempty(rad_idx)
%         continue
%     else
%        for rad_ii = 1:length(rad_idx) 
%         fname_rad = fullfile(strcat(datapath,'rad/',rad_path(rad_idx(rad_ii)).name));
%         temp_rad  = dlmread(fname_rad, ' ',3, 0); 
%         temp_rad_day(rad_ii,:) = temp_rad(:,3);
%        end
%     end
%     
%     time_day = rad_time(rad_idx);
%     
%     spectra = struct('irrad', temp_irrad_day,...
%                      'rad', temp_rad_day,...
%                      'ircoeff', coeff(:,1),...
%                      'rcoeff', coeff(:,2),...
%                      'wl', temp_rad(:,5));
%     WL_range = [745.00,775.00];
%     SIF = SIF_SVD(spectra,WL_range);
%     sif_svd(rad_idx) = SIF.SIF;
%     sif_error(rad_idx) = SIF.SIF_error;
%     sif_day = [time_day,sif_svd(rad_idx)];
%     
% end


% % Store the reflectance
% ref_final       = [time, ref];
% ref_final_time  = sortrows(ref_final,1);
% raw_ref_result  = ref_final_time;
% 
% % Store SIF results
% final_result = [time, f760, rsq, rmse];
% 
% % Time selection
% final_result_time = sortrows(final_result,1);
% 
% raw_final_result = final_result_time;
% 
% %save('SIF760_2014_result.mat');


toc

