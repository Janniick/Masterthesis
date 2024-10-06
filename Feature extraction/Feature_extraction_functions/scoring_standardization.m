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

% Remove rows where any variable has NaN values after conversion
scoring_long = rmmissing(scoring_long);

% % Find the column that contains either 20 or 30 as entries
% % Define the threshold for 90% of the values being either 20 or 30
% threshold = 0.9;
%
% % Check if 90% of the values in 'duration_s' or 'interval_s' are 20 or 30
% if mean(scoring_long.duration_s == 20 | scoring_long.duration_s == 30) >= threshold
%     % If the column with mostly 20/30 values is named 'duration_s'
%     duration_col = 'duration_s';
% elseif mean(scoring_long.interval_s == 20 | scoring_long.interval_s == 30) >= threshold
%     % If the column with mostly 20/30 values is named 'interval_s'
%     duration_col = 'interval_s';
% else
%     error('No column with at least 90% of 20 or 30 durations found.');
% end
%
% % Rename the column to 'duration' and keep only 'epoch', 'stage', and 'duration'
% scoring_standardized = table(scoring_long.epoch, scoring_long.stage, scoring_long.(duration_col), ...
%     'VariableNames', {'epoch', 'stage', 'duration'});
%

%
% % Remove rows where any variable has NaN values after conversion
% scoring_standardized = rmmissing(scoring_standardized);


threshold = 0.9;

% Check if 90% of the values in 'duration_s' or 'interval_s' are 20 or 30
if mean(scoring_long.duration_s == 20 | scoring_long.duration_s == 30) >= threshold
    duration_col = 'duration_s';
elseif mean(scoring_long.interval_s == 20 | scoring_long.interval_s == 30) >= threshold
    duration_col = 'interval_s';
else
    % Merge consecutive 15s
    merge_idx = [];
    duration_col = 'interval_s';
    for i = 1:height(scoring_long) - 1
        if scoring_long.(duration_col)(i) == 15 && scoring_long.(duration_col)(i + 1) == 15
            % Merge the rows
            scoring_long.(duration_col)(i) = 30;  % Update the duration to 30
            merge_idx = [merge_idx; i + 1];  % Mark the second row to be removed
        end
    end

    % Remove the merged rows
    scoring_long(merge_idx, :) = [];

    % Re-enumerate the epochs starting from 1
    scoring_long.epoch = (1:height(scoring_long))';
end

% Find the majority duration (20 or 30 seconds)
majority_duration = mode(scoring_long.(duration_col));

% Find the mismatched rows (where the duration isn't the majority)
mismatched_idx = find(scoring_long.(duration_col) ~= majority_duration);

% Initialize an empty array to hold the new rows
new_rows = [];

% Loop through mismatched indices
i = 1;
while i <= length(mismatched_idx)
    current_idx = mismatched_idx(i);

    % Find consecutive mismatched rows
    j = i;
    while j < length(mismatched_idx) && mismatched_idx(j+1) == mismatched_idx(j) + 1
        j = j + 1;
    end

    % Calculate the total duration of the consecutive mismatched rows
    total_duration = sum(scoring_long.(duration_col)(mismatched_idx(i:j)));

    % Calculate how many new rows we need
    num_new_rows = round(total_duration / majority_duration);

    % Calculate the epoch numbers to replace (incrementing properly)
    replicated_epoch = scoring_long.epoch(mismatched_idx(i)) + (0:(num_new_rows-1))';

    % Replicate the new duration values
    replicated_duration = repmat(majority_duration, num_new_rows, 1);

    % Check if the stages are all the same
    stages = scoring_long.stage(mismatched_idx(i:j));
    if all(stages == stages(1))
        replicated_stage = repelem(stages(1), num_new_rows, 1);  % Same stage
    else
        replicated_stage = linspace(stages(1), stages(end), num_new_rows)';  % Mimic stage progression
    end

    % Create new rows
    new_rows = [new_rows; table(replicated_epoch, replicated_stage, replicated_duration, ...
        'VariableNames', {'epoch', 'stage', 'duration'})];

    % Move to the next set of mismatched rows
    i = j + 1;
end

% Remove the mismatched rows from the original table
scoring_long(mismatched_idx, :) = [];

% Standardize the final table
scoring_standardized = table(scoring_long.epoch, scoring_long.stage, scoring_long.(duration_col), ...
    'VariableNames', {'epoch', 'stage', 'duration'});
if ~isempty(mismatched_idx)

    % Create new rows with the same variables as scoring_long
    new_rows = table(replicated_epoch, replicated_stage, replicated_duration, ...
        'VariableNames', {'epoch', 'stage', 'duration'});  % Use the correct column name f

    % Append the new rows to the original table
    scoring_standardized = [scoring_standardized; new_rows];

    % Sort the table by epoch
    scoring_standardized = sortrows(scoring_standardized, 'epoch');

end

data_type = class(EEG_clean.scoring_long.stage);

% Convert 'stage' to numeric (integers) if they are strings
if data_type == 'string'
    scoring_standardized.stage = str2double(scoring_standardized.stage);  % Convert strings to numeric
end



% Store the standardized scoring table back into EEG_clean
EEG_clean.scoring_long_standardized = scoring_standardized;

% Display the result
disp('Scoring table has been standardized with columns: epoch, stage, and duration.');
end
