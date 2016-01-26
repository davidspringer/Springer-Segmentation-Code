% function states = labelPCGStates(envelope,s1_positions, s2_positions, samplingFrequency, figures)
%
% This function assigns the state labels to a PCG record. 
% This is based on ECG markers, dervied from the R peak and end-T wave locations.
%
%% Inputs:
% envelope: The PCG recording envelope (found in getSchmidtPCGFeatures.m)
% s1_positions: The locations of the R peaks (in samples)
% s2_positions: The locations of the end-T waves (in samples)
% samplingFrequency: The sampling frequency of the PCG recording
% figures (optional): boolean variable dictating the display of figures
%
%% Output:
% states: An array of the state label for each sample in the feature
% vector. The total number of states is 4. Therefore, this is an array of
% values between 1 and 4, such as: [1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,1,1,1],
% illustrating the "true" state label for each sample in the features.
% State 1 = S1 sound
% State 2 = systole
% State 3 = S2 sound
% State 4 = diastole
%
% This code was developed by David Springer for comparison purposes in the
% paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
% where a novel segmentation approach is compared to the paper by Schmidt
% et al:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a 
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
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

function states = labelPCGStates(envelope,s1_positions, s2_positions, samplingFrequency, figures)

if(nargin<5)
    figures = false;
end

states = zeros(length(envelope),1);


%% Timing durations from Schmidt:
mean_S1 = 0.122*samplingFrequency;
std_S1 = 0.022*samplingFrequency;
mean_S2 = 0.092*samplingFrequency;
std_S2 = 0.022*samplingFrequency;

%% Setting the duration from each R-peak to (R-peak+mean_S1) as the first state:
% The R-peak in the ECG coincides with the start of the S1 sound (A. G.
% Tilkian and M. B. Conover, Understanding heart sounds and murmurs: with
% an introduction to lung sounds, 4th ed. Saunders, 2001.)
% Therefore, the duration from each R-peak to the mean_S1 sound duration
% later were labelled as the "true" positions of the S1 sounds:
for i = 1: length(s1_positions)
    %Set an upper bound, incase the window extends over the length of the
    %signal:
    upper_bound = round(min(length(states), s1_positions(i) + mean_S1));
    
    %Set the states between the start of the R peak and the upper bound as
    %state 1:
    states(max([1,s1_positions(i)]):min([upper_bound,length(states)])) = 1;
end

%% Set S2 as state 3 depending on position of end T-wave peak in ECG:
% The second heart sound occurs at approximately the same time as the
% end-T-wave (A. G. Tilkian and M. B. Conover, Understanding heart sounds
% and murmurs: with an introduction to lung sounds, 4th ed. Saunders, 2001.)
% Therefore, for each end-T-wave, find the peak in the envelope around the
% end-T-wave, setting a window centered on this peak as the second heart
% sound state:
for i = 1: length(s2_positions)
    
    %find search window of envelope:
    %T-end +- mean+1sd
    %Set upper and lower bounds, to avoid errors of searching outside size
    %of the signal
    lower_bound = max([s2_positions(i) - floor((mean_S2 + std_S2)),1]);
    upper_bound = min(length(states), ceil(s2_positions(i) + floor(mean_S2 + std_S2)));
    search_window = envelope(lower_bound:upper_bound).*(states(lower_bound:upper_bound)~=1);
    
    % Find the maximum value of the envelope in the search window:
    [~, S2_index] = max(search_window);
    
    %Find the actual index in the envelope of the maximum peak:
    %Make sure this has a max value of the length of the signal:
    S2_index = min(length(states),lower_bound+ S2_index-1);
    
    %Set the states to state 3, centered on the S2 peak, +- 1/2 of the
    %expected S2 sound duration. Again, making sure it does not try to set a
    %value outside of the length of the signal:
    upper_bound = min(length(states), ceil(S2_index +((mean_S2)/2)));
    states(max([ceil(S2_index - ((mean_S2)/2)),1]):upper_bound) = 3;
    
    %Set the spaces between state 3 and the next R peak as state 4:
    if(i<=length(s2_positions))
        %We need to find the next R peak after this S2 sound
        %So, subtract the position of this S2 from the S1 positions
        diffs = (s1_positions - s2_positions(i));
        %Exclude those that are negative (meaning before this S2 occured)
        %by setting them to infinity. They are then excluded when finding
        %the minumum later
        diffs(diffs<0) = inf;
        
        %If the array is empty, then no S1s after this S2, so set to end of
        %signal:
        
        if(isempty(diffs<inf))
            end_pos = length(states);
        else
            %else, send the end position to the minimum diff -1
            [~, index] = min(diffs);
            end_pos = s1_positions(index) -1;
        end
        states(ceil(S2_index +((mean_S2 +(0*std_S2))/2)):end_pos) = 4;
    end
end




%% Setting the first and last sections of the signal
% As all states are derived from either R-peak or end-T-wave locations, the first affirmed state
% in the signal will always be state 1 or state 3. Therefore, until this state, the
% first state should always be set to 4 or 2:

%Find the first step up:
first_location_of_definite_state = find(states ~= 0, 1)-1;

if(first_location_of_definite_state > 1)
    
    if(states(first_location_of_definite_state + 1) == 1)
        states(1:first_location_of_definite_state) = 4;
    end
    
    if(states(first_location_of_definite_state + 1) == 3)
        states(1:first_location_of_definite_state) = 2;
    end
    
end


% Find the last step down:
last_location_of_definite_state = find(states ~= 0, 1,'last');

if(last_location_of_definite_state > 1)
    
    if(states(last_location_of_definite_state) == 1)
        states(last_location_of_definite_state:end) = 2;
    end
    
    if(states(last_location_of_definite_state) == 3)
        states(last_location_of_definite_state:end) = 4;
    end
    
end


states(length(envelope)+1 : end) = [];


%Set everywhere else as state 2:
states(states == 0) = 2;


%% Plotting figures
if(figures)
    figure('Name','Envelope and labelled states');
    plot(envelope);
    hold on;
    plot(states,'r');
    legend('Envelope', 'States');
    pause();
end



