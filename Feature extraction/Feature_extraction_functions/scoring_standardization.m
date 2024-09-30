function EEG_clean = scoring_standardization(EEG_clean)
    % Function to standardize the input scoring table from EEG_clean, 
    % keeping only 'epoch', 'stage', and duration (20 or 30 seconds).
    %
    % Inputs:
    % - EEG_clean: Input structure containing scoring data.
    %
    % Outputs:
    % - EEG_clean: Updated structure with standardized scoring table.

    % Extract the scoring table from EEG_clean
    scoring_long = EEG_clean.scoring_long;
    
    % Find the column that contains either 20 or 30 as entries
    if all(scoring_long.duration_s == 20 | scoring_long.duration_s == 30)
        % If the column with 20/30 is named 'duration_s'
        duration_col = 'duration_s';
    elseif all(scoring_long.interval_s == 20 | scoring_long.interval_s == 30)
        % If the column with 20/30 is named 'interval_s'
        duration_col = 'interval_s';
    else
        error('No valid column with consistent 20 or 30 durations found.');
    end
    
    % Rename the column to 'duration' and keep only 'epoch', 'stage', and 'duration'
    scoring_standardized = table(scoring_long.epoch, scoring_long.stage, scoring_long.(duration_col), ...
        'VariableNames', {'epoch', 'stage', 'duration'});
    
    % Store the standardized scoring table back into EEG_clean
    EEG_clean.scoring_long_standardized = scoring_standardized;
    
    % Display the result
    disp('Scoring table has been standardized with columns: epoch, stage, and duration.');
end
