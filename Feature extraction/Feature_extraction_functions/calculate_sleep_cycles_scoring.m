function [cycle_table, collapsed_values] = calculate_sleep_cycles_scoring(scoring_long_30s)
% Syntax: [cycle_table, collapsed_values] = calculate_sleep_cycles_scoring(scoring_long_30s)
%
% Inputs:
%    Scoring in 30s format
%
% Outputs:
%    cycle_table - MATLAB table containing sleep cycle information:
%                 n_cycle, start_epoch, end_epoch, duration_min.
%    collapsed_values - MATLAB array containing collapsed stages:
%                       [Stage, Start, Length] per row.
% Date:
%    October 2024

%% Step 1: Collapse Consecutive Stages
% Check if scoring_long_30s is empty
if isempty(scoring_long_30s)
    warning('No sleep cycles detected in the hypnogram.');
    cycle_table = table();
    collapsed_values = [];
    return;
end

% Initialize collapsed values storage
collapsed_values = [];
current_stage = scoring_long_30s.stage(1);
current_start = scoring_long_30s.epoch(1);
current_length = scoring_long_30s.duration(1);

% Loop through the rows to collapse consecutive stages
for idx = 2:height(scoring_long_30s)
    % Check if the current stage is the same as the previous one
    if scoring_long_30s.stage(idx) == current_stage
        % Same stage, accumulate the duration
        current_length = current_length + scoring_long_30s.duration(idx);
    else
        % Different stage, store the previous stage details
        collapsed_values = [collapsed_values; current_stage, current_start, current_length / 30];

        % Update to new stage details
        current_stage = scoring_long_30s.stage(idx);
        current_start = scoring_long_30s.epoch(idx);
        current_length = scoring_long_30s.duration(idx);
    end
end

% Add the last stage to collapsed values
collapsed_values = [collapsed_values; current_stage, current_start, current_length / 30];

% Step 1: Filter out rows where duration < 5
collapsed_values_filtered = collapsed_values(collapsed_values(:, 3) >= 6, :);

% Check if there is any remaining data to process
if isempty(collapsed_values_filtered)
    warning('No valid cycles detected after filtering.');
    cycle_table = table();
    collapsed_values = [];
    return;
end

%% Step 2: Re-collapse Consecutive Stages
% Initialize the collapsed values storage for re-collapsed values
recollapsed_values = [];
current_stage = collapsed_values_filtered(1, 1);  % Access the 'stage' column
current_start = collapsed_values_filtered(1, 2);  % Access the 'start' column
current_length = collapsed_values_filtered(1, 3); % Access the 'duration' column

% Loop through the filtered rows to collapse consecutive stages
for idx = 2:size(collapsed_values_filtered, 1)
    if collapsed_values_filtered(idx, 1) == current_stage
        % Same stage, accumulate the duration
        current_length = current_length + collapsed_values_filtered(idx, 3);
    else
        % Different stage, store the previous stage details
        recollapsed_values = [recollapsed_values; current_stage, current_start, current_length];

        % Update to new stage details
        current_stage = collapsed_values_filtered(idx, 1);
        current_start = collapsed_values_filtered(idx, 2);
        current_length = collapsed_values_filtered(idx, 3);
    end
end

% Add the last stage to recollapsed values
recollapsed_values = [recollapsed_values; current_stage, current_start, current_length];

% Update collapsed_values with the new collapsed version
collapsed_values = recollapsed_values;
% Add stop at the end
% Step 1: Remove the duration column from collapsed_values
collapsed_values = collapsed_values(:, 1:2);

% Step 2: Add a final row with 10 and the last epoch index
stop_row = [10, height(scoring_long_30s)];

% Append the final row to collapsed_values
collapsed_values = [collapsed_values; stop_row];

% Initialize the new column
new_column = zeros(size(collapsed_values, 1), 1);

% Populate new column: (second column of next row - 1)
for idx = 1:size(collapsed_values, 1) - 1
    new_column(idx) = collapsed_values(idx + 1, 2) - 1;
end

% For the last row, repeat the value in the second column
new_column(end) = collapsed_values(end, 2);

% Append new column to collapsed_values
collapsed_values = [collapsed_values, new_column];



%% Step 3: Identify Sleep Cycles Based on Stage Transitions
sleep_cycles_struct = struct('start_epoch', {}, 'end_epoch', {}, 'duration_min', {});

start_cycle = 0;
in_rem_cycle = false;  % Flag to mark if in a REM cycle

for idx = 1:size(collapsed_values, 1)
    stage = collapsed_values(idx, 1);  % No conversion needed
    start_epoch = collapsed_values(idx, 2);  % No conversion needed

    if stage == 5  % REM stage
        if start_cycle > 0  % If a NREM start was previously marked
            in_rem_cycle = true;  % Enter REM cycle
        end

    elseif ismember(stage, [1, 2, 3])  % NREM stages
        if in_rem_cycle  % Only end cycle if previously in REM
            % Use the start of the next cycle (NREM) and subtract 1 for end_epoch
            end_cycle = start_epoch - 1;

            % Calculate duration in minutes using start of next epoch minus 1
            duration_min = (end_cycle - start_cycle) * 0.5;

            % Close the current cycle and append to struct
            sleep_cycles_struct(end+1).start_epoch = start_cycle;
            sleep_cycles_struct(end).end_epoch = end_cycle;
            sleep_cycles_struct(end).duration_min = duration_min;

            % Start new cycle
            start_cycle = start_epoch;
            in_rem_cycle = false;  % Reset REM flag for the next cycle
        elseif start_cycle == 0  % Only set start_cycle if it's not set
            start_cycle = start_epoch;
        end
    end
end

% Check for an unfinished cycle at the end of data
if in_rem_cycle && collapsed_values(end, 1) == 10
    % Use the stop epoch minus one as end_epoch
    end_cycle = collapsed_values(end, 2) - 1;

    % Calculate duration in minutes
    duration_min = (end_cycle - start_cycle) * 0.5;

    % Append the final cycle
    sleep_cycles_struct(end+1).start_epoch = start_cycle;
    sleep_cycles_struct(end).end_epoch = end_cycle;
    sleep_cycles_struct(end).duration_min = duration_min;
end

% Convert struct array to table
cycle_table = struct2table(sleep_cycles_struct);

% Add a new column 'n_cycle' with enumerated values
cycle_table.n_cycle = (1:height(cycle_table))';

% Rearrange columns to have 'n_cycle' as the first column
cycle_table = cycle_table(:, [end, 1:end-1]);



end
