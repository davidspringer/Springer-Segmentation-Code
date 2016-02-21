% function [delta, psi, qt] = viterbiDecodePCG_Springer(observation_sequence, pi_vector, b_matrix, total_obs_distribution, heartrate, systolic_time, Fs, figures)
%
% This function calculates the delta, psi and qt matrices associated with
% the Viterbi decoding algorithm from:
% L. R. Rabiner, "A tutorial on hidden Markov models and selected
% applications in speech recognition," Proc. IEEE, vol. 77, no. 2, pp.
% 257-286, Feb. 1989.
% using equations 32a - 35, and equations 68 - 69 to include duration
% dependancy of the states.
%
% This decoding is performed after the observation probabilities have been
% derived from the logistic regression model of Springer et al:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
% Further, this function is extended to allow the duration distributions to extend
% past the beginning and end of the sequence. Without this, the label
% sequence has to start and stop with an "entire" state duration being
% fulfilled. This extension takes away that requirement, by allowing the
% duration distributions to extend past the beginning and end, but only
% considering the observations within the sequence for emission probability
% estimation. More detail can be found in the publication by Springer et
% al., mentioned above.
%
%% Inputs:
% observation_sequence: The observed features
% pi_vector: the array of initial state probabilities, dervived from
% "trainSpringerSegmentationAlgorithm".
% b_matrix: the observation probabilities, dervived from
% "trainSpringerSegmentationAlgorithm".
% heartrate: the heart rate of the PCG, extracted using
% "getHeartRateSchmidt"
% systolic_time: the duration of systole, extracted using
% "getHeartRateSchmidt"
% Fs: the sampling frequency of the observation_sequence
% figures: optional boolean variable to show figures
%
%% Outputs:
% logistic_regression_B_matrix:
% pi_vector:
% total_obs_distribution:
% As Springer et al's algorithm is a duration dependant HMM, there is no
% need to calculate the A_matrix, as the transition between states is only
% dependant on the state durations.
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

function [delta, psi, qt] = viterbiDecodePCG_Springer(observation_sequence, pi_vector, b_matrix, total_obs_distribution, heartrate, systolic_time, Fs,figures)

if nargin < 8
    figures = false;
end

%% Preliminary
springer_options = default_Springer_HSMM_options;

T = length(observation_sequence);
N = 4; % Number of states

% Setting the maximum duration of a single state. This is set to an entire
% heart cycle:
max_duration_D = round((1*(60/heartrate))*Fs);

%Initialising the variables that are needed to find the optimal state path along
%the observation sequence.
%delta_t(j), as defined on page 264 of Rabiner, is the best score (highest
%probability) along a single path, at time t, which accounts for the first
%t observations and ends in State s_j. In this case, the length of the
%matrix is extended by max_duration_D samples, in order to allow the use
%of the extended Viterbi algortithm:
delta = ones(T+ max_duration_D-1,N)*-inf;

%The argument that maximises the transition between states (this is
%basically the previous state that had the highest transition probability
%to the current state) is tracked using the psi variable.
psi = zeros(T+ max_duration_D-1,N);

%An additional variable, that is not included on page 264 or Rabiner, is
%the state duration that maximises the delta variable. This is essential
%for the duration dependant HMM.
psi_duration =zeros(T + max_duration_D-1,N);

%% Setting up observation probs
observation_probs = zeros(T,N);

for n = 1:N
    
    %MLR gives P(state|obs)
    %Therefore, need Bayes to get P(o|state)
    %P(o|state) = P(state|obs) * P(obs) / P(states)
    %Where p(obs) is derived from a MVN distribution from all
    %obserbations, and p(states) is taken from the pi_vector:
    pihat = mnrval(cell2mat(b_matrix(n)),observation_sequence(:,:));
    
    for t = 1:T
        
        Po_correction = mvnpdf(observation_sequence(t,:),cell2mat(total_obs_distribution(1)),cell2mat(total_obs_distribution(2)));
        
        %When saving the coefficients from the logistic
        %regression, it orders them P(class 1) then P(class 2). When
        %training, I label the classes as 0 and 1, so the
        %correct probability would be pihat(2).
        
        observation_probs(t,n) = (pihat(t,2)*Po_correction)/pi_vector(n);
        
    end
end

%% Setting up state duration probabilities, using Gaussian distributions:
[d_distributions, max_S1, min_S1, max_S2, min_S2, max_systole, min_systole, max_diastole, min_diastole] = get_duration_distributions(heartrate,systolic_time);



duration_probs = zeros(N,3*Fs);
duration_sum = zeros(N,1);
for state_j = 1:N
    for d = 1:max_duration_D
        if(state_j == 1)
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_S1 || d > max_S1)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif(state_j==3)
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_S2 || d > max_S2)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif(state_j==2)
            
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_systole|| d > max_systole)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif (state_j==4)
            
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_diastole ||d > max_diastole)
                duration_probs(state_j,d)= realmin;
            end
        end
    end
    duration_sum(state_j) = sum(duration_probs(state_j,:));
end


if(length(duration_probs)>3*Fs)
    duration_probs(:,(3*Fs+1):end) = [];
end

if(figures)
    figure('Name', 'Duration probabilities');
    plot(duration_probs(1,:)./ duration_sum(1),'Linewidth',2);
    hold on;
    plot(duration_probs(2,:)./ duration_sum(2),'r','Linewidth',2);
    hold on;
    plot(duration_probs(3,:)./ duration_sum(3),'g','Linewidth',2);
    hold on;
    plot(duration_probs(4,:)./ duration_sum(4),'k','Linewidth',2);
    hold on;
    legend('S1 Duration','Systolic Duration','S2 Duration','Diastolic Duration');
    pause();
end
%% Perform the actual Viterbi Recursion:


qt = zeros(1,length(delta));
%% Initialisation Step

%Equation 32a and 69a, but leave out the probability of being in
%state i for only 1 sample, as the state could have started before time t =
%0.

delta(1,:) = log(pi_vector) + log(observation_probs(1,:)); %first value is the probability of intially being in each state * probability of observation 1 coming from each state

%Equation 32b
psi(1,:) = -1;


% The state duration probabilities are now used.
%Change the a_matrix to have zeros along the diagonal, therefore, only
%relying on the duration probabilities and observation probabilities to
%influence change in states:
%This would only be valid in sequences where the transition between states
%follows a distinct order.
a_matrix = [0,1,0,0;0 0 1 0; 0 0 0 1;1 0 0 0];


%% Run the core Viterbi algorith

if(springer_options.use_mex)
    
    %% Run Mex code
    % Ensure you have run the mex viterbi_PhysChallenge.c code on the
    % native machine before running this:
    [delta, psi, psi_duration] = viterbi_Springer(N,T,a_matrix,max_duration_D,delta,observation_probs,duration_probs,psi, duration_sum);
    
    
else
    
    %% Recursion
    
    %% The Extended Viterbi algorithm:
    
    %Equations 33a and 33b and 69a, b, c etc:
    %again, ommitting the p(d), as state could have started before t = 1
    
    % This implementation extends the standard implementation of the
    % duration-dependant Viterbi algorithm by allowing the durations to
    % extend beyond the start and end of the time series, thereby allowing
    % states to "start" and "stop" outside of the recorded signal. This
    % addresses the issue of partial states at the beginning and end of the
    % signal being labelled as the incorrect state. For instance, a
    % short-duration diastole at the beginning of a signal looks a lot like
    % systole, and can lead to labelling errors.
    
    % t spans input 2 to T + max_duration_D:
    
    
    for t = 2:T+ max_duration_D-1
        for j = 1:N
            for d = 1:1:max_duration_D
                
                
                %The start of the analysis window, which is the current time
                %step, minus d (the time horizon we are currently looking back),
                %plus 1. The analysis window can be seen to be starting one
                %step back each time the variable d is increased.
                % This is clamped to 1 if extending past the start of the
                % record, and T-1 is extending past the end of the record:
                start_t = t - d;
                if(start_t<1)
                    start_t = 1;
                end
                if(start_t > T-1)
                    start_t = T-1;
                end
                
                %The end of the analysis window, which is the current time
                %step, unless the time has gone past T, the end of the record, in
                %which case it is truncated to T. This allows the analysis
                %window to extend past the end of the record, so that the
                %timing durations of the states do not have to "end" at the end
                %of the record.
                end_t = t;
                if(t>T)
                    end_t = T;
                end
                
                
                %Find the max_delta and index of that from the previous step
                %and the transition to the current step:
                %This is the first half of the expression of equation 33a from
                %Rabiner:
                [max_delta, max_index] = max(delta(start_t,:)+log(a_matrix(:,j))');
                               
                
                %Find the normalised probabilities of the observations over the
                %analysis window:
                probs = prod(observation_probs(start_t:end_t,j));
                
                
                %Find the normalised probabilities of the observations at only
                %the time point at the start of the time window:
                
                if(probs ==0)
                    probs = realmin;
                end
                emission_probs = log(probs);
                
                
                %Keep a running total of the emmission probabilities as the
                %start point of the time window is moved back one step at a
                %time. This is the probability of seeing all the observations
                %in the analysis window in state j:
                
                if(emission_probs == 0 || isnan(emission_probs))
                    emission_probs =realmin;
                end
                
                
                %Find the total probability of transitioning from the last
                %state to this one, with the observations and being in the same
                %state for the analysis window. This is the duration-dependant
                %variation of equation 33a from Rabiner:
                %                 fprintf('log((duration_probs(j,d)./duration_sum(j))):%d\n',log((duration_probs(j,d)./duration_sum(j))));
                delta_temp = max_delta + (emission_probs)+ log((duration_probs(j,d)./duration_sum(j)));
                
                
                %Unlike equation 33a from Rabiner, the maximum delta could come
                %from multiple d values, or from multiple size of the analysis
                %window. Therefore, only keep the maximum delta value over the
                %entire analysis window:
                %If this probability is greater than the last greatest,
                %update the delta matrix and the time duration variable:
                
                
                if(delta_temp>delta(t,j))
                    delta(t,j) = delta_temp;
                    psi(t,j) = max_index;
                    psi_duration(t,j) = d;
                end
                
            end
        end
    end
end


%% Termination

% For the extended case, need to find max prob after end of actual
% sequence:

% Find just the delta after the end of the actual signal
temp_delta = delta(T+1:end,:);
%Find the maximum value in this section, and which state it is in:
[~, idx] = max(temp_delta(:));
[pos, ~] = ind2sub(size(temp_delta), idx);

% Change this position to the real position in delta matrix:
pos = pos+T;

%1) Find the last most probable state
%2) From the psi matrix, find the most likely preceding state
%3) Find the duration of the last state from the psi_duration matrix
%4) From the onset to the offset of this state, set to the most likely state
%5) Repeat steps 2 - 5 until reached the beginning of the signal


%The initial steps 1-4 are equation 34b in Rabiner. 1) finds P*, the most
%likely last state in the sequence, 2) finds the state that precedes the
%last most likely state, 3) finds the onset in time of the last state
%(included due to the duration-dependancy) and 4) sets the most likely last
%state to the q_t variable.

%1)
[~, state] = max(delta(pos,:),[],2);

%2)
offset = pos;
preceding_state = psi(offset,state);

%3)
% state_duration = psi_duration(offset, state);
onset = offset - psi_duration(offset,state)+1;

%4)
qt(onset:offset) = state;

%The state is then updated to the preceding state, found above, which must
%end when the last most likely state started in the observation sequence:
state = preceding_state;

count = 0;
%While the onset of the state is larger than the maximum duration
%specified:
while(onset > 2)
    
    %2)
    offset = onset-1;
    %     offset_array(offset,1) = inf;
    preceding_state = psi(offset,state);
    %     offset_array(offset,2) = preceding_state;
    
    
    %3)
    %     state_duration = psi_duration(offset, state);
    onset = offset - psi_duration(offset,state)+1;
    
    %4)
    %     offset_array(onset:offset,3) = state;
    
    if(onset<2)
        onset = 1;
    end
    qt(onset:offset) = state;
    state = preceding_state;
    count = count +1;
    
    if(count> 1000)
        break;
    end
end

qt = qt(1:T);


