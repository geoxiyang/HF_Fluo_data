%% TEST CALCULATION OF SIF

% BF_WF_file = dir('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr15_Brown_Test/BareFiber_WhiteReference_*');
% BF_PL_file = dir('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr15_Brown_Test/BareFiber_Plants_*');
% 
% [calib,a,b] = importdata('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOfile/Apr152015_Calibration_CosineCorrector_OOIIrrad.cal_OOIIrrad.cal','\t',9);
% [calib2,a2,b2] = importdata('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOfile/Apr072015_Calibration_OOIIrrad.cal_OOIIrrad.cal','\t',9);
% 
% 
% for ii = 1: 50
%     
%     
%     REF_F  = fullfile('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr15_Brown_Test/',BF_WF_file(ii).name);
%     CAL_F  = fullfile('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr15_Brown_Test/',BF_PL_file(ii).name);
%     
%     [ref,del,header]    = importdata(REF_F,'\t',16);
%     [cal,del2,header2]  = importdata(CAL_F,'\t',16);
% 
%     % Conversion for uncalibrated file
% %     xdata = calib.data(:,2).*ref.data(:,2)*10/((6e-1*3.9*3.9*1e-2*0.5*(calib.data(2,1)-calib.data(1,1))));
% %     ydata = calib2.data(:,2).*cal.data(:,2)*10/((1.50000E-1*1.0*1.0*1e-2*0.5*(calib2.data(2,1)-calib2.data(1,1)))*(pi()*sind(12.5)^2));
%     
%     % Conversion for calibrated file
%     xdata = ref.data(:,2)*10/(pi()*sind(12.5)^2);
%     ydata = cal.data(:,2)*10/(pi()*sind(12.5)^2);
% 
% plotyy(ref.data(:,1),xdata,ref.data(:,1),ydata);
%     
%     WL_range = [745.00,780.00]; % Including O2A: 717.00 780.00; 745.00 780.00
%     irrad = xdata';
%     rad   = ydata';
%     
%     coeff = ones(1044,2);
%     
%     
%     spectra = struct('irrad', irrad,...
%                      'rad', rad,...
%                      'ircoeff', coeff(:,1),...
%                      'rcoeff', coeff(:,2),...
%                      'wl', ref.data(:,1));
%     WL_range = [745.00,780.00]; % Including O2A: 717.00 780.00; 745.00 780.00
%     SIF_result = SIF_SVD(spectra,WL_range,2,1,0.02,6);
%     sif_svd(ii) = SIF_result.SIF;
%     sif_error(ii) = SIF_result.SIF_error;
%     sif_rel(ii) = SIF_result.SIF_relative;
% 
% end
% 
% plot(1:1:length(sif_svd),sif_rel,'ko');
% opath = '/Users/xiyang/Desktop/barefilber_plants.png';
% print(gcf, '-dpng','-r300', opath);

cd('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr29_Brown_test')

    REF_F  = fullfile('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr29_Brown_Test/QEPR000001_10-45-34-328.txt');
    PAL_F  = fullfile('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOdata/Apr29_Brown_Test/QEPR000001_10-45-42-877.txt');
    
    Ref_int=1.5e-1;
    Pal_int=1.5e-1;
    Cal_int=1.65e-1;
    
    
    [ref,del,header]    = importdata(REF_F,'\t',16);
    [pal,del2,header2]  = importdata(PAL_F,'\t',16);
    [calib2,a2,b2]      = importdata('/Volumes/XiYangResearch/Data/13.Fluorescence/2015test/OOfile/Apr072015_Calibration_OOIIrrad.cal_OOIIrrad.cal','\t',9);

    % Conversion for uncalibrated file
    xdata = calib2.data(:,2).*ref.data(:,2)*10/((Ref_int*1.0*1.0*1e-2*0.5*(calib2.data(2,1)-calib2.data(1,1)))*(pi()*sind(12.5)^2));
    ydata = calib2.data(:,2).*pal.data(:,2)*10/((Pal_int*1.0*1.0*1e-2*0.5*(calib2.data(2,1)-calib2.data(1,1)))*(pi()*sind(12.5)^2));
    
   plot(ref.data(:,1),xdata(:),'r-')
   xlim([730,780])
   hold on
   plot(pal.data(:,1),ydata(:),'b-')
   hold off
    
   % plotyy(ref.data(:,1),xdata(:),pal.data(:,1),ydata(:));
    %xlim([740,780]);
    WL_range = [745.00,780.00]; % Including O2A: 717.00 780.00; 745.00 780.00
    irrad = xdata';
    rad   = ydata';
    
    coeff = ones(1044,2);
    
    
    spectra = struct('irrad', irrad,...
                     'rad', rad,...
                     'ircoeff', coeff(:,1),...
                     'rcoeff', coeff(:,2),...
                     'wl', ref.data(:,1));
    WL_range = [745.00,780.00]; % Including O2A: 717.00 780.00; 745.00 780.00
    SIF_result = SIF_SVD(spectra,WL_range,2,1,0.02,6);
    sif_svd = SIF_result.SIF;
    sif_error = SIF_result.SIF_error;
    sif_rel = SIF_result.SIF_relative;
