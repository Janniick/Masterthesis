%% Function to Expand Condensed Scoring for Different Formats
function EEG_clean = scoring_expansion(EEG_clean)

    % Detect which type of scoring format is present
    scoring = EEG_clean.scoring;  % Assuming EEG_clean.scoring is the table
    
    % Check for the presence of different column headers to detect the format
    if all(ismember({'epoch', 'stage', 'duration_s', 'sample_rate_hz'}, scoring.Properties.VariableNames))
        % ACJ format
        fprintf('Detected ACJ format.\n');
        
        % Initialize an empty table for the expanded epochs
        expanded_scoring = table([], [], [], [], 'VariableNames', scoring.Properties.VariableNames);

        % Loop through each row of the original scoring table
        for i = 1:height(scoring)
            current_epoch = scoring.epoch(i);
            if i < height(scoring)
                next_epoch = scoring.epoch(i + 1);
            else
                next_epoch = current_epoch;  % Handle the last row
            end
            % Expand all epochs between current and next_epoch
            for epoch_num = current_epoch:(next_epoch - 1)
                new_row = table(epoch_num, scoring.stage(i), scoring.duration_s(i), scoring.sample_rate_hz(i), ...
                    'VariableNames', scoring.Properties.VariableNames);
                expanded_scoring = [expanded_scoring; new_row];  % Append to the expanded table
            end
        end

    elseif all(ismember({'nr', 'segment', 'epoch', 'starttime_s', 'duration_s', 'interval_s', 'stage'}, scoring.Properties.VariableNames))
        % ASD format
        fprintf('Detected ASD format.\n');

        % Initialize an empty table for the expanded epochs
        expanded_scoring = table([], [], [], [], [], [], [], 'VariableNames', scoring.Properties.VariableNames);

        % Loop through each row of the original scoring table
        for i = 1:height(scoring)
            current_epoch = scoring.epoch(i);
            if i < height(scoring)
                next_epoch = scoring.epoch(i + 1);
            else
                next_epoch = current_epoch;  % Handle the last row
            end
            % Expand all epochs between current and next_epoch
            for epoch_num = current_epoch:(next_epoch - 1)
                new_row = table(scoring.nr(i), scoring.segment(i), epoch_num, scoring.starttime_s(i), ...
                    scoring.duration_s(i), scoring.interval_s(i), scoring.stage(i), ...
                    'VariableNames', scoring.Properties.VariableNames);
                expanded_scoring = [expanded_scoring; new_row];  % Append to the expanded table
            end
        end

    else
        error('Unknown scoring format.');
    end

    % Assign the expanded table back to EEG_clean.scoring_long
    EEG_clean.scoring_long = expanded_scoring;
    fprintf('Scoring expansion completed.\n');
end
