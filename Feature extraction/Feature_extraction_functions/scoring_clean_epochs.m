function EEG_clean = scoring_clean_epochs(EEG_clean)
    % Function to extract rows corresponding to clean epochs from EEG_clean.scoring_long_30s.
    % An epoch is considered clean if it contains no 1's in EEG_clean.rem_ind_all data.
    %
    % Inputs:
    % - EEG_clean: EEG data structure with fields 'scoring_long_30s', 'rem_ind_all', and 'srate'.
    %
    % Outputs:
    % - EEG_clean: EEG data structure with the 'clean_epochs' field containing only clean epochs from scoring_long_30s.
    
    % Number of samples per epoch
    samples_per_epoch = EEG_clean.scoring_long_30s.duration(1) * EEG_clean.srate;
    
    % Number of epochs in EEG_clean.scoring_long_30s
    num_epochs = height(EEG_clean.scoring_long_30s);
    
    % Initialize array to store clean epochs
    clean_epoch_info = [];
    
    % Loop through each epoch and check if it contains any bad signal (1's)
    for i = 1:num_epochs
        % Get the start and end indices for this epoch in EEG_clean.rem_ind_all
        start_idx = EEG_clean.scoring_long_30s.start_index(i);
        end_idx = EEG_clean.scoring_long_30s.end_index(i);
        
        % Extract the data for this epoch from rem_ind_all
        epoch_data = EEG_clean.rem_ind_all(start_idx:end_idx);
        
        % If there are no 1's in this epoch, store the corresponding row from scoring_long_30s
        if all(epoch_data == 0)
            clean_epoch_info = [clean_epoch_info; EEG_clean.scoring_long_30s(i, :)];
        end
    end
    
    % Assign the rows of clean epochs back to EEG_clean as a new field
    EEG_clean.clean_epochs = clean_epoch_info;
    
    % Calculate the percentage of clean epochs
    clean_epochs_count = size(clean_epoch_info, 1);
    clean_epochs_percent = (clean_epochs_count / num_epochs) * 100;
    
    % Display the number and percentage of clean epochs
    disp(['Number of clean epochs: ', num2str(clean_epochs_count)]);
    disp(['Percentage of clean epochs: ', num2str(clean_epochs_percent), '%']);
end
