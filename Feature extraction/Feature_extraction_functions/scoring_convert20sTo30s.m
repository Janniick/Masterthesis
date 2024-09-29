% Function: scoring_convert20sTo30s
%
% Description:
% This function converts a sleep scoring table with 20-second epochs into a new table with 30-second epochs. 
% If the input table is already in 30-second epochs, no changes are made. The function works by combining every 
% three consecutive 20-second epochs into two 30-second epochs, keeping the stage data from the first and third 
% epochs, respectively.
%
% Additionally, the function calculates the start and end sample indices for each 30-second epoch based on the 
% sample rate and adds them as new columns to the output table.
%
% Inputs:
% - scoring_long: A table containing sleep scoring data with columns for epoch number, sleep stage, 
%   duration (in seconds), sample rate, and optionally start and end sample indices.
%
% Outputs:
% - scoring_long_30s: A new table with 30-second epochs, updated epoch numbers, and calculated start and end 
%   sample indices for each 30-second epoch.
%
% Usage Example:
% scoring_long_30s = scoring_convert20sTo30s(scoring_long)
% 
% This function assumes that all the epochs in the input table are either 20 or 30 seconds long and will return an 
% error if the input table contains inconsistent durations. If the epochs
% are already in 30s format, the function won't do anything
%
% Author: Jannick Mauron
% Date: 2024
%


function scoring_long_30s = scoring_convert20sTo30s(scoring_long)
    % Check if the file is already in 30-second epochs
    if all(scoring_long.duration_s == 30)
        disp('The file is already in 30-second epochs. No changes made.');
        scoring_long_30s = scoring_long;  % Return the original table
        return;
    end
    
    % Check if the epochs are in 20 seconds format
    if all(scoring_long.duration_s == 20)
        % Initialize a new table for the 30s epochs
        scoring_long_30s = table([], [], [], [], 'VariableNames', scoring_long.Properties.VariableNames);
        
        % Loop through every set of 3 consecutive 20s epochs
        for i = 1:3:height(scoring_long)
            % Check if there are at least 3 consecutive epochs remaining
            if i + 2 > height(scoring_long)
                % If fewer than 3 epochs remain, discard them
                break;
            end
            
            % Create the first new 30s epoch by keeping the first stage
            new_row_1 = table(scoring_long.epoch(i), scoring_long.stage(i), 30, scoring_long.sample_rate_hz(i), ...
                'VariableNames', scoring_long.Properties.VariableNames);
            
            % Create the second new 30s epoch by keeping the third stage
            new_row_2 = table(scoring_long.epoch(i+1), scoring_long.stage(i+2), 30, scoring_long.sample_rate_hz(i), ...
                'VariableNames', scoring_long.Properties.VariableNames);
            
            % Append the new rows to the 30s scoring table
            scoring_long_30s = [scoring_long_30s; new_row_1; new_row_2];
            

        % After all epochs have been processed, renumber the epochs consecutively
        num_epochs = height(scoring_long_30s);
        scoring_long_30s.epoch = (1:num_epochs)';  % Renumber epochs from 1 to x
        end
    else
        error('The file contains epochs with inconsistent durations.');
    end
    
   % Add start_index and end_index of each respective epoch
    samples_per_epoch = scoring_long_30s.duration_s(1) * scoring_long_30s.sample_rate_hz(1);  % Number of samples per 30s epoch
    start_indices = (0:num_epochs-1)' * samples_per_epoch + 1;  % Start indices for each epoch
    end_indices = start_indices + samples_per_epoch - 1;  % End indices for each epoch

    % Add the start_index and end_index as new columns to the scoring_long_30s table
    scoring_long_30s.start_index = start_indices;
    scoring_long_30s.end_index = end_indices;

    % Display the new 30s scoring table
    disp('Converted to 30-second epochs');
end