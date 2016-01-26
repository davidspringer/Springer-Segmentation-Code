% function low_pass_filtered_signal = butterworth_low_pass_filter(original_signal,order,cutoff,sampling_frequency, figures)
%
% Low-pass filter a given signal using a forward-backward, zero-phase
% butterworth low-pass filter.
%
%% INPUTS:
% original_signal: The 1D signal to be filtered
% order: The order of the filter (1,2,3,4 etc). NOTE: This order is
% effectively doubled as this function uses a forward-backward filter that
% ensures zero phase distortion
% cutoff: The frequency cutoff for the low-pass filter (in Hz)
% sampling_frequency: The sampling frequency of the signal being filtered
% (in Hz).
% figures (optional): boolean variable dictating the display of figures
%
%% OUTPUTS:
% low_pass_filtered_signal: the low-pass filtered signal.
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., ?Logistic Regression-HSMM-based Heart Sound
% Segmentation,? IEEE Trans. Biomed. Eng., In Press, 2015.
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

function low_pass_filtered_signal = butterworth_low_pass_filter(original_signal,order,cutoff,sampling_frequency, figures)

if nargin < 5,
    figures = 0;
end

%Get the butterworth filter coefficients
[B_low,A_low] = butter(order,2*cutoff/sampling_frequency,'low');

if(figures)
    figure('Name','Low-pass filter frequency response');
    [sos,g] = zp2sos(B_low,A_low,1);	     % Convert to SOS form
    Hd = dfilt.df2tsos(sos,g);   % Create a dfilt object
    h = fvtool(Hd);	             % Plot magnitude response
    set(h,'Analysis','freq')	     % Display frequency response
end


%Forward-backward filter the original signal using the butterworth
%coefficients, ensuring zero phase distortion
low_pass_filtered_signal = filtfilt(B_low,A_low,original_signal);

if(figures)
    figure('Name','Original vs. low-pass filtered signal');
    plot(original_signal);
    hold on;
    plot(low_pass_filtered_signal,'r');
    legend('Original Signal', 'Low-pass filtered signal');
    pause();
end