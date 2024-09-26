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
    
    % Display the new 30s scoring table
    disp('Converted to 30-second epochs');
end