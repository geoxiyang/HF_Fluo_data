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

%  Determine index range based on WL_range
fl_index = 


































end