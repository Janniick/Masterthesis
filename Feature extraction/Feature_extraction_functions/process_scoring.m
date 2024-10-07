function EEG_clean = process_scoring(EEG_clean)
    % Expand the scoring if it's condensed
    % Generates EEG_clean.scoring_long
    EEG_clean = scoring_expansion(EEG_clean);
    % Standardize the scoring
    % Only keep epoch, stage and duration
    % Generates EEG_clean.scoring_long_standardized
    EEG_clean = scoring_standardization(EEG_clean);
    % convert 20s to 30s duration if needed
    % Generates EEG_clean.scoring_long_30s
    EEG_clean = scoring_convert20sTo30s(EEG_clean);
    % add indexes of every epoch to the scoring
    % Updates EEG_clean.scoring_long_30s
    EEG_clean = scoring_add_start_end_indices(EEG_clean);
end
