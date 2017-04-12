function SIF = SIF_SVD(spectra, WL_range, polyDegree1, polyDegree2, variance_threshold, min_rank, numFLv, var_m, var_b)
%% Introduction
% Calculate SIF from spectra collected from field spectrometers.
% Revised for my analysis of Harvard Forest data from Ari Kornfeld
% (Carnegie Institute for Science)'s code.
% Also the SVD part is from Christian Frankenberg and Luis Guanter
% History: 
% v1.0    11/19/2014   Xi Yang (xi_yang@brown.edu)
% Reference:
% Guanter, L., M. Rossini, R. Colombo, M. Meroni, C. Frankenberg, J.-E. Lee
% , and J. Joiner (2013), Using field spectroscopy to assess the potential 
% of statistical approaches for the retrieval of sun-induced chlorophyll 
% fluorescence from ground and space, Remote Sensing of Environment, 133, 
% doi:10.1016/j.rse.2013.01.017.

%% Usage
%  1. Input
%       spectra             -- a structure in which
%         spectra.irrad     -- irradiance DN
%         spectra.rad       -- radiance DN
%         spectra.ircoeff   -- irradiance coefficient
%         spectra.rcoeff    -- radiance coefficient
%         spectra.wl        -- wavelength
%       WL_range            -- wavelength range of the region to extract
%                           SIF signal;                
%       polyDegree1         --
%       polyDegree2         --
%       variance_threshold  --
%       min_rank            --
%       numFLv              --
%       var_m               --
%       var_b               --
%  2. Output
%       SIF                 --
%  3. Example
%       SIF = SIF_SVD()
%

%% Start of the program %%
%                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Setting the defaults

if nargin < 3
    polyDegree1 = 2;
end
if nargin < 4
    polyDegree2 = 1;
end
if nargin < 5
    variance_threshold = 1;
end
if nargin < 6
    min_rank = 2;
end

if nargin < 7
    numFLv = 1;
end

if nargin < 8
    var_m = 0; var_b = 1;
end


%% 2. Prepare spectra for SIF

%  Determine index range based on WL_range, fl_index == TRUE if it is in
%  that range
fl_index = selectRanges(spectra.wl,WL_range);

%  Prepare the fluorescence shape array
%  This dataset is measured by Ari Kornfeld using several plants in
%  Stanford campus
%  FLraw: 1st column is wavelength, 2:8 are SVs (?)
FLraw       = load('/Volumes/XiYangBackUp/src/SIFviewer v2 source/FLshape3.mat');
%  FL_SV: interpolate the spectrum of fluorescence shape into current
%  wavelength range
FL_SV       = interp1(FLraw.FLshape3(:,1),FLraw.FLshape3(:,2:8),spectra.wl','pchip');
% Test
% irrCoefs    = spectra.ircoeff;
%  Convert fluorescence shape (in radiance units) into DN (should I use irradiance or radiance coeff?)
% Test
% FL_spec     = MxV(@rdivide, FL_SV(:, 1:numFLv), irrCoefs) ;
%  ???????
FL760idx    = selectRanges(spectra.wl, [759.5, 760.5]);
% the values of the SVs at 760 nm (rescaled to the spectrometer sensitivity)
% Test
% FL760       = mean(FL_spec(FL760idx, :), 1)'; 
% FL760coef   = mean(irrCoefs(FL760idx) );
FL760       = mean(FL_SV(FL760idx, 1:numFLv), 1)';

%  Prepare the input (setup) that is used in the training and fitting of
%  SVD method
%  pixels: the indices of wavelength that is used in the fitting;
%  wl: wavelength (nm);
%  variance_threshold: The threshold of the choice of Eigenvalues
%  polynomial degree for broadband features (for Eigenvector 1)
%  polynomial degree for broadband features (for Eigenvector 2)
setup.pixels            = fl_index;
setup.wl                = spectra.wl(fl_index);
setup.variance_threshold=variance_threshold; %default: 1    (0.15); 
setup.polyDegree1       = polyDegree1; % default 2?
setup.polyDegree2       = polyDegree2; % default 1
setup.hf                = FL_SV(fl_index, numFLv);
% Test
% setup.irrCoefs          = irrCoefs(fl_index);
% setup.hf                = FL_spec(fl_index, numFLv);


%% 3. Calculate SIF
irrad = MxV(@times,spectra(:).irrad', spectra.ircoeff);
rad   = MxV(@times,spectra(:).rad',spectra.rcoeff);
irrad = irrad';
rad   = rad';

% Test
% irrad = spectra(:).irrad;
% rad   = spectra(:).rad;

%  Compute SVD from irradiance spectra (DN)
TrainingResult = trainSVD_SIF(irrad(:,fl_index),setup,min_rank);

%  Compute SVD from radiance spectra
%  First select wavelength (for reference?)
WL758 = selectRanges(spectra.wl, [757.65 758.06]);
%  scale FL760 to match the scaling of K (??)
FL760 = FL760./TrainingResult.Ksd( TrainingResult.FLidx)';
%  now calculate SIF
SIF   = solveSVD_SIF(rad, TrainingResult, FL760, WL758, var_m, var_b);



%% End of the program   %%
%                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%


end
function R = trainSVD_SIF(spec, setup, min_rank)
%trainSVD_SIF   Computes singular value composition of a set of spectra and computes the resulting fit matrix needed for subsequent Solar Induced Fluorescence fitting.
%   trainSVD_SIF(spec, setup)
%   Input Arguments:
%     spec - the spectrometer data corresponding to setup.wl
%     setup:
%       setup.wavlengths: wavelengths to use
%       setup.variance_threshold: variance threshold (in %) of eigenvalues,
%          (i.e. determines amount of Eigenvectors to be used)
%       setup.polyDegree1: degree of polynomial used for EV1
%       setup.polyDegree2: degree of polynomial used for EV1
%     min_rank: the minimum # of eigenvectors to use (1 or 2)
%
%   Christian Frankenberg, with small mods by Ari Kornfeld (JAK)
%
%   $Revision: 1.0 $  $Date: 2013/02/07 $
%    JAK Mod   2013-05-01
%        Transpose V & R.wl to match revised shape of wavelength vector (it's now a row instead of column)
%   JAK Mod   2013-05-29 rearrange setup to handle arbitrary pixel specification

% get final wavelength range:
%R.pixels = setup.pixels;  % pixel mask, can be used to subset the full dataset later
%R.wl = setup.wavelengths';
R = setup;  % preserve the setup data
R.meanwl = mean(R.wl);
% Perform SVD of the dataset in spec
[U,S,V] = svd(spec, 'econ');
% Normalize Eigenvalues of SVD to compute fraction of explained variance
% and use variance_threshold to determine how many Eigenvectors to include
% in the analysis
R.norm_S = diag(S)/sum(diag(S));
ind = find(R.norm_S*100>setup.variance_threshold);
if(length(ind)<min_rank)
    % provide at least min_rank Eigenvectors (typ. 2) as long as they are actually available
    ind = 1:min(min_rank, size(V, 2));
end

% Return at least the first 20 Eigenvectors passing the eigenvalue threshold (for display in SIFviewer):
evToreturn = 1:min( max(20, max(ind)), size(V, 2)); % at least 10, but only if there are at least 10 EV! 
R.EV = -V(:,evToreturn);  % Eigenvectors always come out -ve, so reverse them here.
R.U = -U(:, evToreturn);
R.S = S(evToreturn, evToreturn);
R.norm_S = R.norm_S(evToreturn);
% R.Snorm = R.norm_S;

% Compute Matrices for the Eigenvectors convolved with polynomial terms in
% lambda space (change it here a little, make it as function of
% (lambda-mean(lambda))
%  JAK: reorganize to match the order used by polyval()
idx = 1;
for i=setup.polyDegree1:-1:0
    R.K_poly1(:,idx) =  R.EV(:, 1).*(R.wl - R.meanwl).^i;  
    idx = idx + 1;
end  

R.K_poly2=[];
if ( length(ind) > 1 )
    idx = 1;
    for i=setup.polyDegree2:-1:0
        R.K_poly2(:,idx) =  R.EV(: , 2).*(R.wl - R.meanwl).^i; 
        idx = idx + 1;
    end
else
    setup.polyDegree2 = -1; % (note degree 0 is a valid input. Also need -1 for R.poly2idx, below)
end

if (max(ind) > 2 )
    EV3plus =  R.EV(: , 3:max(ind));
else
    EV3plus = [];
end
% SIF shape vector, hf, setup in the calling function, above (JAK)

% set up final Jacobian matrix used for linear inversion:
%  (JAK note: zero-order terms are now included in the polynomials, so only need those beyond two)
R.K = [R.K_poly1, R.K_poly2, EV3plus, setup.hf];
% return indeces of the polynomials in the Jacobian matrix
R.poly1idx = 1:size(R.K_poly1, 2); % note: polynomial includes the zero-order term
R.poly2idx = R.poly1idx(end) + (1:size(R.K_poly2, 2)); % can be empty
R.nEV = max(ind); % total # EV in the Jacobian
R.otherEVidx = max([R.poly1idx, R.poly2idx]) + (1:size(EV3plus, 2)); % can be empty
R.FLidx = (size(R.K, 2) - size(setup.hf, 2)+1):size(R.K, 2);

R.Ksd = std(R.K); % to disable, use: ones(1, size(R.K, 2)); %
%R.K = R.K * diag(sparse( 1./R.Ksd )); % normalize each column so SD = 1, to improve rank, least-square, calculations
R.K = MxV(@rdivide, R.K, R.Ksd ); % normalize each column so SD = 1, to improve rank, least-square, calculations

end % function trainSVD_SIF

function R = solveSVD_SIF(full_spec, train, FL760, WL758, var_m, var_b)  %, pixels
%solveSFD_SIF   Computing solar induced fluorescence (SIF) using the SVD approach.
%   solveSVD_SIF(spec, K, pixels, var_m, var_b)
%   performs the least squares fit of all spectra in "spec" (just using the
%   index range provided in "pixels" and using "K" as Jacobi Matrix
%   (previously obtained by using a training dataset). The method largely
%   follows Luis' paper in RSE.
%  Optional: var_m and var_b allow weighting the least squares by the measured values.
%  in the Avantes spectrometer, for example, variance = 0.65 * y + 3600
%   where y is the measured "counts";  Default: var_m = 0; var_b = 1
%
%   Returns structure R:
%   R.SIF (SIF)
%   R.SIF_relative (fractional SIF related to continuum level radiance; in %)
%   R.sv  (entire state vector of the fit (dimesnion depends on provided K,
%   for simplicity, last element should be SIF)
%   R.mod  (modeled spectrum)
%   R.meas (measured spectrum)
%   R.rms  (fit residuum RMS)
%   R.var  (variance of fit residual)
%   R.SNR  (estimated SNR of the scene; based on the fit RMS)
%   R.SIF_error  (estimated 1-sigma error of the SIF fit)
%
%   Christian Frankenberg
%
%   $Revision: 1.0 $  $Date: 2013/02/07 $

if nargin < 5
    var_m = 0; var_b = 1;
end

K = train.K;
spec = full_spec(:, train.pixels);

% for all spectra at once:
    if var_m == 0
    % perform linear unweighted least squares
        result = K\spec';
    else
    % sadly, weighting by Y cannot be done in one line since weights differ for each y[i]:
        result = zeros(size(K, 2), size(spec, 1)); % initialize result
        for idx = 1:size(spec, 1)
            % weighted least squares:  X'WX Beta = X'Wy 
            % use "shortcut" X" = wX, y" = wy where w = sqrt(W)  
            w = sqrt(spec(idx, :).*var_m + var_b); % variance is a linear function of signal strength
            result(:, idx) = (diag(1./w) *K)\ (spec(idx, :) ./ w)';
            % For some reason the following results in ill-conditioned matrices, whereas the
            %  code above does not.
            %w = (spec(idx, :).*var_m + var_b); % variance is a linear function of signal strength
            %W1 = diag(max(w)./w);
            %result(:, idx) = (K' *W1 *K) \ (K' * W1 *spec(idx, :)')  ;
        end
    end
    % save all state vector elements as well
    R.sv = result';
    % Save SIF result (allow for multiple fluorescent shape vectors
    numFLv = length(FL760); % the input values of the flourescence eigenvectors at 760
    %  the x values used for regression & FL760 have been scaled by their standard deviation, 
    R.SIF = R.sv(:, (end-numFLv+1):end) * FL760 ;
    % relative SIF the original led to confusing results. Instead, always normalize to a fixed wavelength
    R.SIF_relative = 100 * R.SIF ./ mean(full_spec(:,WL758), 2);  % make relative to the same spot regardless of SIF frame [was: max(spec(i,end)) ]
    %R.SIF_relative = 100 * R.SIF ./ spec(:,selectRanges(spec.wavelengths, [745 746]));  % was: max(spec(i,end)) - not max(spec(i, :))
    % Compute RMS of the fit:
    % modeled spectrum:
    R.mod = (K*result)';
    % measured spectrum:
    R.meas = spec;
    % fit RMS:
    R.rms = rms(R.meas - R.mod, 2);  % rms(spec(i,:)'-K*result);
    % fit variance:
    R.var = var(R.meas - R.mod, 0, 2);   %var(spec(i,:)'-K*result);
    % Compute SNR from data itself:
    R.SNR = mean(spec, 2)./R.rms;
    % Compute fit error (variance) based on noise estimated from residuals:
    % estimate error on CH4 scaling factor:
    %Se = R.var(i)*inv(transpose(K)*K);
    %R.SIF_error(i) = sqrt(Se(end,end)); 
    Se1 = inv(transpose(K)*K); % label? K'K is the covariance matrix;  K'y (K'K)^1 is least-squares fit
    R.SIF_error = sqrt(R.var * Se1(end,end));

end

function notrun()

for i=1:length(spec(:,1))
    % perform linear (unweighted for now) least squares
    res = K\spec(i,:)';
    % Save methane result 
    R.SIF(i) = res(end);
    R.SIF_relative(i) = res(end)/max(spec(i,end))*100;
    % save all state vector elements as well
    R.sv(i,:) = res;
    % Compute RMS of the fit:
    % modeled spectrum:
    R.mod(i,:) = K*res;
    % measured spectrum:
    R.meas(i,:) = spec(i,:);
    % fit RMS:
    R.rms(i) = rms(spec(i,:)'-K*res);
    % fit variance:
    R.var(i) = var(spec(i,:)'-K*res);
    % Compute SNR from data itself:
    R.SNR(i) = mean(spec(i,:))/R.rms(i);
    % Compute fit error (variance) based on noise estimated from residuals:
    Se = R.var(i)*inv(transpose(K)*K);
    % estimate error on CH4 scaling factor:
    R.SIF_error(i) = sqrt(Se(end,end)); 
    
end

end %function solveSVD_SIF

function F = fspectrum(wavelengths) 
if false
     % not very good out of range
flpoly = [1.84932187469354e-22     -6.23627727252016e-21 ...
     -1.47836464392334e-18      4.55784730072008e-17 ...
      4.55896682723519e-15     -1.17847092764609e-13 ...
      -6.9789723765354e-12      1.23909167959185e-10 ...
      5.57778419606746e-09     -3.29795084331999e-08 ...
     -1.48325449617467e-06     -3.81635447631006e-05 ...
      -0.00119647445081813        0.0177309326271857 ...
          2.32894181199785];
F = polyval(flpoly, wavelengths-730);
else
    % Coefficients (with 95% confidence bounds):
       a1 =        2.14; %  (2.121, 2.158)
       b1 =       682.6; %  (682.5, 682.6)
       c1 =        11.2; %  (11.08, 11.32)
       a2 =       60.32; %  (-5.325e+06, 5.325e+06)
       b2 =       737.8; %  (427.6, 1048)
       c2 =       11.36; %  (-949.2, 971.9)
       a3 =      -59.89; %  (-5.325e+06, 5.325e+06)
       b3 =       737.8; %  (426.2, 1049)
       c3 =       11.33; %  (-954.8, 977.4)
       a4 =        2.03; %  (2.009, 2.05)
       b4 =       727.3; %  (727, 727.6)
       c4 =       48.98; %  (48.63, 49.34)
       x = wavelengths;
     F = a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2) + ... 
              a3*exp(-((x-b3)/c3).^2) + a4*exp(-((x-b4)/c4).^2);
end

end

































