% function springer_options = default_Springer_HSMM_options()
%
% The default options to be used with the Springer segmentation algorithm.
% USAGE: springer_options = default_Springer_HSMM_options
%
% Developed for use in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function springer_options = default_Springer_HSMM_options()

%% The sampling frequency at which to extract signal features:
springer_options.audio_Fs = 1000;

%% The downsampled frequency
%Set to 50 in Springer paper
springer_options.audio_segmentation_Fs = 50;


%% Tolerance for S1 and S2 localization
springer_options.segmentation_tolerance = 0.1;%seconds

%% Whether to use the mex code or not:
% The mex code currently has a bug. This will be fixed asap.
springer_options.use_mex = false;

%% Whether to use the wavelet function or not:
springer_options.include_wavelet_feature = false;

