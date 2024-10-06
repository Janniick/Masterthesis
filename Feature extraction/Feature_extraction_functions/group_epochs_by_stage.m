function stage_epochs = group_epochs_by_stage(clean_epochs)
    % group_epochs_by_stage - Group epochs by sleep stages.
    %
    % Syntax: stage_epochs = group_epochs_by_stage(clean_epochs)
    %
    % Inputs:
    %    clean_epochs - Struct containing epoch data with fields:
    %                  epoch, stage, start_index, end_index.
    %
    % Outputs:
    %    stage_epochs - Struct containing grouped start and end indices:
    %                   NREM1, NREM2, NREM3, REM.
    %
    % Example:
    %    stage_epochs = group_epochs_by_stage(EEG_clean.clean_epochs);
    %
    % Author:
    %    Jannick (Adapted by ChatGPT)
    %
    % Date:
    %    October 2024

    % Extract necessary data
    epochs = clean_epochs.epoch;
    stages = clean_epochs.stage;
    start_indices = clean_epochs.start_index;
    end_indices = clean_epochs.end_index;
    
    % Initialize a struct to hold start and end indices grouped by stages
    stage_epochs = struct('NREM1', [], 'NREM2', [], 'NREM3', [], 'REM', []);
    
    % Loop through all epochs and group start and end indices by their stage
    for i = 1:length(stages)
        % Ignore stage 0 (wake) and other undefined stages
        if stages(i) == 1
            stage_epochs.NREM1 = [stage_epochs.NREM1; start_indices(i), end_indices(i)];
        elseif stages(i) == 2
            stage_epochs.NREM2 = [stage_epochs.NREM2; start_indices(i), end_indices(i)];
        elseif stages(i) == 3 || stages(i) == 4
            stage_epochs.NREM3 = [stage_epochs.NREM3; start_indices(i), end_indices(i)];
        elseif stages(i) == 5
            stage_epochs.REM = [stage_epochs.REM; start_indices(i), end_indices(i)];
        end
    end
end
