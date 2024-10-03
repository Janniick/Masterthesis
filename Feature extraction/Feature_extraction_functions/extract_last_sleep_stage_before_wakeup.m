function last_sleep_stage = extract_last_sleep_stage_before_wakeup(hypnogram)
    % Define the number of epochs to check (e.g., 6)
    epochs_to_check = 6;
    
    % Find the last occurrence of a switch from a non-zero stage to stage 0
    wake_transition_indices = find(diff(hypnogram == 0) == 1);
    
    if isempty(wake_transition_indices)
        % If no wake transition is found, use the last 'epochs_to_check' epochs
        last_epochs = hypnogram(max(1, end-epochs_to_check+1):end);
    else
        % Get the index of the last transition to wake (stage 0)
        last_wake_transition = wake_transition_indices(end);
        
        % Extract the last 'epochs_to_check' epochs before the last transition to wake
        last_epochs = hypnogram(max(1, last_wake_transition-epochs_to_check+1):last_wake_transition);
    end
    
    % Find the sleep stages (excluding 0)
    sleep_stages = last_epochs(last_epochs ~= 0);
    
    if isempty(sleep_stages)
        last_sleep_stage = NaN;  % No non-zero stages found
    else
        % Find the frequency of each stage
        [unique_stages, ~, stage_indices] = unique(sleep_stages);
        stage_counts = histc(stage_indices, 1:numel(unique_stages));
        
        % Check if there is a clear winner (a single most frequent stage)
        max_count = max(stage_counts);
        if sum(stage_counts == max_count) > 1
            % If there's a tie, take the last occurring stage
            last_sleep_stage = sleep_stages(end);
        else
            % Otherwise, take the most frequent stage
            last_sleep_stage = unique_stages(stage_counts == max_count);
        end
    end
end

