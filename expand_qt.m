% function expanded_qt = expand_qt(original_qt, old_fs, new_fs, new_length)
% 
% Function to expand the derived HMM states to a higher sampling frequency. 
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% original_qt: the original derived states from the HMM
% old_fs: the old sampling frequency of the original_qt
% new_fs: the desired sampling frequency
% new_length: the desired length of the qt signal

%% Outputs:
% expanded_qt: the expanded qt, to the new FS and length
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

function expanded_qt = expand_qt(original_qt, old_fs, new_fs, new_length)

original_qt = original_qt(:)';
expanded_qt = zeros(new_length,1);

indeces_of_changes = find(diff(original_qt));

indeces_of_changes = [indeces_of_changes, length(original_qt)];

start_index = 0;
for i = 1:length(indeces_of_changes)
    
    start_index;
    end_index = indeces_of_changes(i);
    
    mid_point = round((end_index - start_index)/2) + start_index;
    
    value_at_mid_point = original_qt(mid_point);
    
    expanded_start_index = round((start_index./old_fs).*new_fs) + 1;
    expanded_end_index = round((end_index./(old_fs)).*new_fs);
    
    if(expanded_end_index > new_length)
        expanded_end_index = new_length;
    end
    
    expanded_qt(expanded_start_index:expanded_end_index) = value_at_mid_point;

    start_index = end_index;
end