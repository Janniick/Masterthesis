function EEG_clean = scoring_add_start_end_indices(EEG_clean)
    % Function to add start_index and end_index for each epoch in EEG_clean.scoring_long_30s.
    %
    % Inputs:
    % - EEG_clean: A struct that contains the scoring table (EEG_clean.scoring_long_30s) 
    %              and the sample rate (EEG_clean.srate).
    %
    % Outputs:
    % - EEG_clean: The input struct with the updated scoring table, 
    %              now including 'start_index' and 'end_index' columns.
    
    % Extract the scoring table and sample rate
    scoring_table = EEG_clean.scoring_long_30s;
    sample_rate = EEG_clean.srate;

    % Number of epochs
    num_epochs = height(scoring_table);
    
    % Calculate samples per epoch
    samples_per_epoch = scoring_table.duration(1) * sample_rate;  % Number of samples per epoch
    
    % Calculate start and end indices
    start_indices = (0:num_epochs-1)' * samples_per_epoch + 1;  % Start indices for each epoch
    end_indices = start_indices + samples_per_epoch - 1;        % End indices for each epoch

    % Add the start_index and end_index as new columns to the scoring table
    scoring_table.start_index = start_indices;
    scoring_table.end_index = end_indices;

    % Update the EEG_clean struct with the modified scoring table
    EEG_clean.scoring_long_30s = scoring_table;

    % Display a message
    disp('Added start_index and end_index to EEG_clean.scoring_long_30s.');
end
