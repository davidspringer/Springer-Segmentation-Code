% function [heartRate systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)
%
% Derive the heart rate and the sytolic time interval from a PCG recording.
% This is used in the duration-dependant HMM-based segmentation of the PCG
% recording.
%
% This method is based on analysis of the autocorrelation function, and the
% positions of the peaks therein.
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a 
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% audio_data: The raw audio data from the PCG recording
% Fs: the sampling frequency of the audio recording
% figures: optional boolean to display figures
%
%% OUTPUTS:
% heartRate: the heart rate of the PCG in beats per minute
% systolicTimeInterval: the duration of systole, as derived from the
% autocorrelation function, in seconds
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

function [heartRate, systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)

if nargin < 3
    figures = false;
end

%% Get heatrate:
% From Schmidt:
% "The duration of the heart cycle is estimated as the time from lag zero
% to the highest peaks between 500 and 2000 ms in the resulting
% autocorrelation"
% This is performed after filtering and spike removal:

%% 25-400Hz 4th order Butterworth band pass
audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs, false);
audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs);

%% Spike removal from the original paper:
audio_data = schmidt_spike_removal(audio_data,Fs);

%% Find the homomorphic envelope
homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);

%% Find the autocorrelation:

signal_autocorrelation = autocorr(homomorphic_envelope,length(homomorphic_envelope)-1, [] , 2);

min_index = 0.5*Fs;
max_index = 2*Fs;


[~, index] = max(signal_autocorrelation(min_index:max_index));
true_index = index+min_index-1;

heartRate = 60/(true_index/Fs);


%% Find the systolic time interval:
% From Schmidt: "The systolic duration is defined as the time from lag zero
% to the highest peak in the interval between 200 ms and half of the heart
% cycle duration"


max_sys_duration = round(((60/heartRate)*Fs)/2);
min_sys_duration = round(0.2*Fs);

[~, pos] = max(signal_autocorrelation(min_sys_duration:max_sys_duration));
systolicTimeInterval = (min_sys_duration+pos)/Fs;


if(figures)
    figure('Name', 'Heart rate calculation figure');
    plot(signal_autocorrelation);
    hold on;
    plot(true_index, signal_autocorrelation(true_index),'ro');
    plot((min_sys_duration+pos), signal_autocorrelation((min_sys_duration+pos)), 'mo');
    xlabel('Samples');
    legend('Autocorrelation', 'Position of max peak used to calculate HR', 'Position of max peak within systolic interval');
end


