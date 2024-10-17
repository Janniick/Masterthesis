function ECG_results = calculate_ECG_by_cycles_and_stages(EEG_clean, cycle_table, collapsed_values)
    % Initialize the results structure
    ECG_results = struct();
    
    % Channels to process (ECG data)
    nr_cycles = height(cycle_table);  % Number of sleep cycles
    nr_epochs = height(EEG_clean.scoring_long_30s);  % Total number of epochs
    
    % Initialize storage for NREM and REM data for each cycle
    first_cycle_data = struct();
    last_cycle_data = struct();
    
    % Initialize storage for all cycles
    all_NREM_data = [];
    all_REM_data = [];
    
    % Loop through each sleep cycle
    for cycle_idx = 1:nr_cycles
        cycle_name = ['Cycle_' num2str(cycle_idx)];
        start_epoch = cycle_table.start_epoch(cycle_idx);
        end_epoch = cycle_table.end_epoch(cycle_idx);
        
        % Initialize structures for NREM and REM results for this cycle
        NREM_data = [];
        REM_data = [];
        
        % Loop through the collapsed values and check sleep stages within this cycle
        for i = 1:size(collapsed_values, 1)
            stage = collapsed_values(i, 1);  % Sleep stage
            epoch_start = collapsed_values(i, 2);  % Start epoch of this stage

            % Only process if the stage is within the current cycle range
            if epoch_start >= start_epoch && epoch_start <= end_epoch
                epoch_end = epoch_start + collapsed_values(i, 3) - 1;  % End epoch
                epoch_indices = epoch_start:epoch_end;

                if stage == 1 || stage == 2 || stage == 3  % NREM stages (1 = NREM1, 2 = NREM2, 3 = NREM3)
                    % Calculate ECG for NREM stages
                    for idx = epoch_indices
                        if idx <= nr_epochs
                            indexs = EEG_clean.scoring_long_30s{idx, 'start_index'}:EEG_clean.scoring_long_30s{idx, 'end_index'};
                            ECG_values = calculate_ECG_values(EEG_clean, indexs);
                            NREM_data = [NREM_data; struct2array(ECG_values)];  % Collect ECG data for NREM
                        end
                    end
                elseif stage == 5  % REM stage
                    % Calculate ECG for REM stages
                    for idx = epoch_indices
                        if idx <= nr_epochs
                            indexs = EEG_clean.scoring_long_30s{idx, 'start_index'}:EEG_clean.scoring_long_30s{idx, 'end_index'};
                            ECG_values = calculate_ECG_values(EEG_clean, indexs);
                            REM_data = [REM_data; struct2array(ECG_values)];  % Collect ECG data for REM
                        end
                    end
                end
            end
        end

        % Collect all cycle data for overall mean calculations
        all_NREM_data = [all_NREM_data; NREM_data];
        all_REM_data = [all_REM_data; REM_data];

        % Calculate mean and standard deviation for NREM and REM
        if ~isempty(NREM_data)
            NREM_mean = mean(NREM_data, 1, 'omitnan');
            NREM_std = std(NREM_data, 0, 1, 'omitnan');
            ECG_results.(cycle_name).NREM.mean = struct('ecg_mean', NREM_mean(1), 'ecg_rmssd', NREM_mean(2), ...
                'ecg_sdnn', NREM_mean(3), 'ecg_sdsd', NREM_mean(4), 'ecg_pnn50', NREM_mean(5), ...
                'ecg_LF', NREM_mean(6), 'ecg_HF', NREM_mean(7), 'ecg_att_en', NREM_mean(8), 'br_permin', NREM_mean(9));
            ECG_results.(cycle_name).NREM.std = struct('ecg_mean', NREM_std(1), 'ecg_LF', NREM_std(6), ...
                'ecg_HF', NREM_std(7), 'ecg_att_en', NREM_std(8), 'br_permin', NREM_std(9));
        end
        
        if ~isempty(REM_data)
            REM_mean = mean(REM_data, 1, 'omitnan');
            REM_std = std(REM_data, 0, 1, 'omitnan');
            ECG_results.(cycle_name).REM.mean = struct('ecg_mean', REM_mean(1), 'ecg_rmssd', REM_mean(2), ...
                'ecg_sdnn', REM_mean(3), 'ecg_sdsd', REM_mean(4), 'ecg_pnn50', REM_mean(5), ...
                'ecg_LF', REM_mean(6), 'ecg_HF', REM_mean(7), 'ecg_att_en', REM_mean(8), 'br_permin', REM_mean(9));
            ECG_results.(cycle_name).REM.std = struct('ecg_mean', REM_std(1), 'ecg_LF', REM_std(6), ...
                'ecg_HF', REM_std(7), 'ecg_att_en', REM_std(8), 'br_permin', REM_std(9));
        end

        % Calculate differences between NREM and REM
        if ~isempty(NREM_data) && ~isempty(REM_data)
            ECG_results.(cycle_name).NREM_REM_differences = struct( ...
                'ecg_mean_diff', ECG_results.(cycle_name).NREM.mean.ecg_mean - ECG_results.(cycle_name).REM.mean.ecg_mean, ...
                'ecg_att_en_diff', ECG_results.(cycle_name).NREM.mean.ecg_att_en - ECG_results.(cycle_name).REM.mean.ecg_att_en, ...
                'br_permin_diff', ECG_results.(cycle_name).NREM.mean.br_permin - ECG_results.(cycle_name).REM.mean.br_permin);
        end

        % Store the first and last cycle data for later comparison
        if cycle_idx == 1
            first_cycle_data.NREM_mean = ECG_results.(cycle_name).NREM.mean;
            first_cycle_data.REM_mean = ECG_results.(cycle_name).REM.mean;
        elseif cycle_idx == nr_cycles
            last_cycle_data.NREM_mean = ECG_results.(cycle_name).NREM.mean;
            last_cycle_data.REM_mean = ECG_results.(cycle_name).REM.mean;
        end
    end

    % Calculate differences between first and last cycle for each measure
    if ~isempty(first_cycle_data) && ~isempty(last_cycle_data)
        ECG_results.Cycle_First_Last_Differences = struct( ...
            'NREM_ecg_mean_diff', last_cycle_data.NREM_mean.ecg_mean - first_cycle_data.NREM_mean.ecg_mean, ...
            'NREM_ecg_att_en_diff', last_cycle_data.NREM_mean.ecg_att_en - first_cycle_data.NREM_mean.ecg_att_en, ...
            'NREM_br_permin_diff', last_cycle_data.NREM_mean.br_permin - first_cycle_data.NREM_mean.br_permin, ...
            'REM_ecg_mean_diff', last_cycle_data.REM_mean.ecg_mean - first_cycle_data.REM_mean.ecg_mean, ...
            'REM_ecg_att_en_diff', last_cycle_data.REM_mean.ecg_att_en - first_cycle_data.REM_mean.ecg_att_en, ...
            'REM_br_permin_diff', last_cycle_data.REM_mean.br_permin - first_cycle_data.REM_mean.br_permin);
    end

    % Calculate overall means and standard deviations across all cycles
    if ~isempty(all_NREM_data)
        all_NREM_mean = mean(all_NREM_data, 1, 'omitnan');
        all_NREM_std = std(all_NREM_data, 0, 1, 'omitnan');
        ECG_results.Overall.NREM.mean = struct('ecg_mean', all_NREM_mean(1), 'ecg_rmssd', all_NREM_mean(2), ...
            'ecg_sdnn', all_NREM_mean(3), 'ecg_sdsd', all_NREM_mean(4), 'ecg_pnn50', all_NREM_mean(5), ...
            'ecg_LF', all_NREM_mean(6), 'ecg_HF', all_NREM_mean(7), 'ecg_att_en', all_NREM_mean(8), 'br_permin', all_NREM_mean(9));
        ECG_results.Overall.NREM.std = struct('ecg_mean', all_NREM_std(1), 'ecg_LF', all_NREM_std(6), ...
            'ecg_HF', all_NREM_std(7), 'ecg_att_en', all_NREM_std(8), 'br_permin', all_NREM_std(9));
    end

    if ~isempty(all_REM_data)
        all_REM_mean = mean(all_REM_data, 1, 'omitnan');
        all_REM_std = std(all_REM_data, 0, 1, 'omitnan');
        ECG_results.Overall.REM.mean = struct('ecg_mean', all_REM_mean(1), 'ecg_rmssd', all_REM_mean(2), ...
            'ecg_sdnn', all_REM_mean(3), 'ecg_sdsd', all_REM_mean(4), 'ecg_pnn50', all_REM_mean(5), ...
            'ecg_LF', all_REM_mean(6), 'ecg_HF', all_REM_mean(7), 'ecg_att_en', all_REM_mean(8), 'br_permin', all_REM_mean(9));
        ECG_results.Overall.REM.std = struct('ecg_mean', all_REM_std(1), 'ecg_LF', all_REM_std(6), ...
            'ecg_HF', all_REM_std(7), 'ecg_att_en', all_REM_std(8), 'br_permin', all_REM_std(9));
        ECG_results.Overall.NREM_REM_differences = struct('ecg_mean_diff', ECG_results.Overall.NREM.mean.ecg_mean - ECG_results.Overall.REM.mean.ecg_mean, ...
                'ecg_att_en_diff', ECG_results.Overall.NREM.mean.ecg_att_en - ECG_results.Overall.REM.mean.ecg_att_en, ...
                'br_permin_diff', ECG_results.Overall.NREM.mean.br_permin - ECG_results.Overall.REM.mean.br_permin);
    end
end
