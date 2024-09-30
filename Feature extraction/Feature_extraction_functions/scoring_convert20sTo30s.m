function EEG_clean = convert20sto30s(EEG_clean)
    % Function to convert a standardized scoring table in EEG_clean from
    % 20-second epochs to 30-second epochs or leave it unchanged if it is already in 30s.
    %
    % Inputs:
    % - EEG_clean: Input structure containing the standardized scoring table.
    %
    % Outputs:
    % - EEG_clean: Updated structure with the new 30-second epoch table.

    % Extract the standardized scoring table
    scoring_standardized = EEG_clean.scoring_long_standardized;
    
    % Check if the epochs are already in 30-second format
    if all(scoring_standardized.duration == 30)
        disp('The file is already in 30-second epochs. No changes made.');
        EEG_clean.scoring_long_30s = scoring_standardized;
        return;
    elseif all(scoring_standardized.duration == 20)
        % Initialize a new table for the 30s epochs
        scoring_30s = table([], [], [], 'VariableNames', {'epoch', 'stage', 'duration'});
        
        % Loop through every set of 3 consecutive 20s epochs
        for i = 1:3:height(scoring_standardized)
            % Check if there are at least 3 consecutive epochs remaining
            if i + 2 > height(scoring_standardized)
                % If fewer than 3 epochs remain, discard them
                break;
            end
            
            % Create the first new 30s epoch by keeping the first stage
            new_row_1 = table(scoring_standardized.epoch(i), scoring_standardized.stage(i), 30, ...
                'VariableNames', {'epoch', 'stage', 'duration'});
            
            % Create the second new 30s epoch by keeping the third stage
            new_row_2 = table(scoring_standardized.epoch(i+1), scoring_standardized.stage(i+2), 30, ...
                'VariableNames', {'epoch', 'stage', 'duration'});
            
            % Append the new rows to the 30s scoring table
            scoring_30s = [scoring_30s; new_row_1; new_row_2];
        end
        
        % After all epochs have been processed, renumber the epochs consecutively
        num_epochs = height(scoring_30s);
        scoring_30s.epoch = (1:num_epochs)';  % Renumber epochs from 1 to x
        
        % Store the 30s scoring table in EEG_clean
        EEG_clean.scoring_long_30s = scoring_30s;
        
    else
        error('The file contains epochs with inconsistent durations.');
    end
    
    % Display the result
    disp('Converted to 30-second epochs.');
end
