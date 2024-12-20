%% add functions
addpath(genpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Feature extraction\Feature_extraction_functions'));
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
py.importlib.import_module('yasa');

%% load python
pyenv('Version', 'C:\Users\janni\AppData\Local\Programs\Python\Python311\python.EXE');
% pyrun("import yasa", "yasa");
%% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data';  % Replace with your folder path
% Search for all files ending with '_preprocessed_EOG.mat' in all subdirectories
all_mat_files = dir(fullfile(base_folder, '**', '*_preprocessed_EOG.mat'));
%% Load existing results_table if it exists
results_table_file = fullfile(base_folder, 'results_table.mat');  % Save in base_folder
if exist(results_table_file, 'file')
    load(results_table_file, 'results_table');
else
    results_table = table();
end

% Load existing error_log if it exists
error_log_file = fullfile(base_folder, 'error_log.mat');  % Save in base_folder
if exist(error_log_file, 'file')
    load(error_log_file, 'error_log');
else
    error_log = table('Size', [0, 4], ...
                  'VariableTypes', {'double', 'string', 'string', 'string'}, ...
                  'VariableNames', {'i', 'study', 'participant', 'error_message'});
end

%% bands
delta_b = [0.5, 4];
theta_b = [4, 8];
alpha_b = [8, 12];
beta_b = [15, 30];
gamma_b = [30, 40];
spindle_b = [12, 16];

for idx = 1:length(all_mat_files)
    try
        tic;
        % extract the name of the file
        filename = all_mat_files(idx).name;
        disp(['------------------------ Analysing File ', num2str(idx), ' of ', num2str(length(all_mat_files)), ': ', filename]);
        % Split the filename by underscores ('_') and periods ('.')
        split_parts = split(filename, {'_', '.'});
        % Extract the study and participant variables
        study = split_parts{1};
        participant = split_parts{2};
        % Check if the study is 'MOG'
        if strcmp(study, 'MOG')
            % Append split_parts{3} to participant
            participant = strcat(participant, '_', split_parts{3});
        end

        % Check if the file has already been processed
        if ~isempty(results_table)
            existing_entries = strcmp(results_table.study, study) & strcmp(results_table.participant, participant);
            if any(existing_entries)
                disp(['File ', filename, ' already processed. Skipping.']);
                continue; % Skip to next file
            end
        end

        % Check if the file has already encountered an error in error_log
        if exist('error_log', 'var') && ~isempty(error_log)
            error_entries = strcmp(error_log.study, study) & strcmp(error_log.participant, participant);
            if any(error_entries)
                disp(['File ', filename, ' already logged with an error. Skipping.']);
                continue; % Skip to next file
            end
        end
        
        % load the file
        load([all_mat_files(idx).folder, '\', all_mat_files(idx).name]);

        % Initialize results_row as a struct
        results_row = struct();

        % Save study and participant
        results_row.i = idx;
        results_row.study = study;
        results_row.participant = participant;


        % initialize things for spectrum and so on
        srate = EEG_clean.srate;                 % signal rate
        duration_window = 4*srate;               % 4 second windows
        which_window = hanning(duration_window); % hanning
        overlap_window = duration_window*0.5;    % 50%
        % multitaper spectrogram
        frequency_range = [0 40]; % Limit frequencies from 0 to 40 Hz
        window_s = 5; % window length in seconds
        step_s = 2; % step size of moving the window along
        freq_res = 1; % resolution in HZ, higher the number the smoother
        tw = window_s*(1/2); % window length* (frequency resolution)
        nr_tapers = floor(2*(tw))-1;
        taper_params = [tw nr_tapers]; % Time half-bandwidth and number of tapers
        window_params = [window_s step_s]; % Window size s with step size of 2s
        min_nfft = 0; % No minimum nfft
        detrend_opt = 'off'; % detrend each window by subtracting the average
        weighting = 'unity'; % weight each taper at 1
        plot_on = false; % plot spectrogram
        verbose = false; % print extra info



        % Make the scoring uniform
        EEG_clean = process_scoring(EEG_clean);



        % Calculate data and scoring durations and cut the longer one accordingly
        % Assuming that they both have the same starting time point
        data_duration_samples = length(EEG_clean.data);  % Total samples in data
        scoring_duration_samples = height(EEG_clean.scoring_long_30s) * 30 * EEG_clean.srate;  % Scoring duration in samples
        % Display the duration in hours
        data_duration_hours = data_duration_samples / EEG_clean.srate / 3600;
        scoring_duration_hours = scoring_duration_samples / EEG_clean.srate / 3600;
        disp(['Duration of EEG_clean.data (hours): ', num2str(data_duration_hours)]);
        disp(['Duration of EEG_clean.scoring_long_30s (hours): ', num2str(scoring_duration_hours)]);
        % Check if data or scoring is longer and cut the longer one
        if data_duration_samples > scoring_duration_samples
            % Data is longer, so cut it
            EEG_clean.data = EEG_clean.data(:, 1:scoring_duration_samples);
            EEG_clean.times = EEG_clean.times(:, 1:scoring_duration_samples);
            EEG_clean.data_ECG = EEG_clean.data_ECG(:, 1:scoring_duration_samples);
            EEG_clean.data_EOG = EEG_clean.data_EOG(:, 1:scoring_duration_samples);
            EEG_clean.rem_ind_all = EEG_clean.rem_ind_all(:, 1:scoring_duration_samples);
            disp(['Data length reduced to match scoring duration: ', num2str(scoring_duration_hours), ' hours']);
        elseif scoring_duration_samples > data_duration_samples
            % Scoring is longer, so cut it
            scoring_rows = floor(data_duration_samples / EEG_clean.srate / 30);  % Number of 30s epochs to match data
            EEG_clean.scoring_long_30s = EEG_clean.scoring_long_30s(1:scoring_rows, :);
            disp(['Scoring length reduced to match data duration: ', num2str(data_duration_hours), ' hours']);
        else
            disp('Data and scoring durations already match.');
        end



        % Making the hypnogram for the whole night
        % Extract the sleep stages from the whole night
        sleep_stages_whole = EEG_clean.scoring_long_30s.stage;
        hypnogram = create_hypnogram(sleep_stages_whole);
        % Convert MATLAB hypnogram to Python format using py.numpy.array
        hypnogram_py = py.numpy.array(hypnogram);
        % Upsample the whole night hypnogram
        epoch_duration_sec = 30;
        repeats_per_epoch = epoch_duration_sec * EEG_clean.srate;  % Number of samples per epoch
        hypnogram_upsampled = repelem(hypnogram, repeats_per_epoch);  % Upsample the hypnogram
        hypnogram_upsampled_py = py.numpy.array(hypnogram_upsampled);  % Convert to Python numpy array



        % Calculate sleep cycles
        [cycle_table, collapsed_values] = calculate_sleep_cycles_scoring(EEG_clean.scoring_long_30s);
        % Save the amount of cycles
        results_row.n_cycles = cycle_table.n_cycle(end);


        % Check for clean epochs by checking if no disturbance was found in rem_ind_all
        % Generates EEG_clean.clean_epochs
        EEG_clean = scoring_clean_epochs(EEG_clean);



        % Extract only clean epoch data based on EEG_clean.clean_epochs
        % Initialize an empty array to store clean data
        clean_data = [];
        % Loop through each clean epoch and concatenate the corresponding data
        for i = 1:height(EEG_clean.clean_epochs)
            start_idx = EEG_clean.clean_epochs.start_index(i);  % Get start index for current epoch
            end_idx = EEG_clean.clean_epochs.end_index(i);      % Get end index for current epoch
            clean_data = [clean_data, EEG_clean.data(:, start_idx:end_idx)];  % Concatenate clean data
        end
        % Store clean data in a new field in EEG_clean
        EEG_clean.clean_data = clean_data;
        % Do the same for times
        clean_times = [];
        % Loop through each clean epoch and concatenate the corresponding data
        for i = 1:height(EEG_clean.clean_epochs)
            start_idx = EEG_clean.clean_epochs.start_index(i);  % Get start index for current epoch
            end_idx = EEG_clean.clean_epochs.end_index(i);      % Get end index for current epoch
            clean_times = [clean_times, EEG_clean.times(:, start_idx:end_idx)];  % Concatenate clean data
        end
        % Store clean data in a new field in EEG_clean
        EEG_clean.clean_times = clean_times;
        % Convert the clean EEG data to Python-friendly format (NumPy array)
        clean_data_for_py = cell(1, size(EEG_clean.clean_data, 1));
        for row = 1:size(EEG_clean.clean_data, 1)
            clean_data_for_py{row} = EEG_clean.clean_data(row, :);
        end
        clean_data_for_py = py.numpy.array(clean_data_for_py);




        % Prepare hypnogram data of clean epochs
        sleep_stages_clean = EEG_clean.clean_epochs.stage;
        % Determine the number of epochs (each epoch is 30 seconds)
        num_epochs = length(sleep_stages_clean);
        % Initialize the hypnogram vector to store sleep stages per epoch
        hypnogram_clean = zeros(1, num_epochs);
        % Loop over each epoch and assign the corresponding sleep stage to the hypnogram
        % 0=wake, 1=NREM1, 2=NREM2, 3&4=NREM3, 5=REM, 6<=movement
        for epoch = 1:num_epochs
            stage = sleep_stages_clean(epoch);
            if stage == 3 || stage == 4
                % Combine stage 3 and 4 into NREM3 (using 3 as NREM3 code)
                hypnogram_clean(epoch) = 3;
            elseif stage >= 6
                % Anything above 5 (movement) is set to 6
                hypnogram_clean(epoch) = 6;
            else
                % Use the stage as it is for other sleep stages (0 to 5)
                hypnogram_clean(epoch) = stage;
            end
        end
        % Upsample hypnogram_clean
        epoch_duration_sec = 30;  % Duration of each epoch in seconds
        srate = EEG_clean.srate;  % Sampling rate of the EEG data
        % Repeat each entry in hypnogram_clean for 30 * srate times
        repeats_per_epoch = epoch_duration_sec * srate;  % Number of samples per epoch
        hypnogram_clean_upsampled = repelem(hypnogram_clean, repeats_per_epoch);


        %% Extraction of whole night features
        % Call YASA's sleep statistics function
        disp('- Extract whole night sleep statistics');
        % Slightly modify hypnogram for yasa function
        new_hypnogram = hypnogram;
        new_hypnogram(new_hypnogram == 5) = 4;
        new_hypnogram_py = py.numpy.array(new_hypnogram);
        sleep_stats = py.yasa.sleep_statistics(new_hypnogram_py, sf_hyp=1/30);
        sleep_stats = dictionary(sleep_stats);
        % Save the generated values
        sleep_stat_keys = keys(sleep_stats);
        for i = 1:length(sleep_stat_keys)
            % Get the key
            key = sleep_stat_keys{i};  % Original key from the dictionary
            % Replace '%' with 'perc_' and ensure the name is valid
            if startsWith(key, '%')
                valid_key = ['perc_' key(2:end)];  % Replace '%' with 'perc_' only at the start
            else
                valid_key = key;
            end
            valid_key = matlab.lang.makeValidName(valid_key);  % Ensure the modified key is valid as a field name
            % Get the corresponding value
            value = sleep_stats(key);
            % Assign to the structure with the modified key as the field name
            results_row.(valid_key) = value;
        end



        % Calculate stage transitions
        stage_transitions = py.yasa.transition_matrix(hypnogram_py);
        % Extract counts and probabilities from the result
        counts = stage_transitions{1};  % Transition counts
        probs = stage_transitions{2};   % Transition probabilities
        % Convert counts and probs to NumPy arrays
        counts_np = counts.to_numpy();
        probs_np = probs.to_numpy();
        % Convert Python counts and probabilities to MATLAB arrays
        counts_mat = double(counts_np);  % Convert to MATLAB numeric matrix
        probs_mat = double(probs_np);    % Convert to MATLAB numeric matrix
        % Define the stage numbers
        stages = [0, 1, 2, 3, 5, 6];
        % Loop through the transition probabilities matrix
        for row = 1:size(probs_mat, 1)
            for col = 1:size(probs_mat, 2)
                % Construct a valid field name using the specific stage numbers
                field_name = sprintf('prob_%d_to_%d', stages(row), stages(col));
                % Store the transition probability in results_row
                results_row.(field_name) = probs_mat(row, col);
            end
        end
        % Loop through the transition count matrix to store the values
        for row = 1:size(counts_mat, 1)
            for col = 1:size(counts_mat, 2)
                % Construct a valid field name for the transition
                field_name = sprintf('count_%d_to_%d', stages(row), stages(col));
                % Store the transition probability in results_row
                results_row.(field_name) = counts_mat(row, col);
            end
        end



        % Extract the last sleep stage before wakeup
        last_stage_before_wakeup = extract_last_sleep_stage_before_wakeup(hypnogram);
        % Save
        results_row.last_stage_before_wakeup = last_stage_before_wakeup;



        % wasserstein interhemispheric spectral measure
        % This function computes the Wasserstein distance between two EEG channels based on their power spectral density (PSD).
        % The PSD is calculated using the Welch method over the specified window and overlap parameters. The function then calculates
        % the cumulative distribution of the PSD in the frequency range of 0.5 to 40 Hz and compares the two channels using the
        % Wasserstein distance metric.
        % done over whole night
        % remove no channels
        disp('- Calculate wasserstein distance');
        remove_channel_which = zeros(1,8);
        duration_window = 4*EEG_clean.srate;
        which_window = hanning(duration_window);
        overlap_window = duration_window*0.5;
        % Saves it automatically with appropriate name, e.g. wasserstein_f3f4
        wasserstein_distance(EEG_clean, 'F3', 'F4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'C3', 'C4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'P3', 'P4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'O1', 'O2', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        results_row.wasserstein_F3F4 = wasserstein_f3f4;
        results_row.wasserstein_C3C4 = wasserstein_c3c4;
        results_row.wasserstein_P3P4 = wasserstein_p3p4;
        results_row.wasserstein_O1O2 = wasserstein_o1o2;



        % Calculate ECG values during sleep
        disp('- Calculate ECG values');
        ECG_results = calculate_ECG_by_cycles_and_stages(EEG_clean, cycle_table, collapsed_values);
        % Save
        results_row = addNestedFieldsToResultsRow(ECG_results.Cycle_First_Last_Differences, 'ECG_results.Cycle_First_Last_Differences', results_row);
        results_row = addNestedFieldsToResultsRow(ECG_results.Overall, 'ECG.Overall', results_row);

        

        % Calculate complexity measures
        complexity_measures = extract_complexity_measures(EEG_clean, srate, delta_b, theta_b, alpha_b, beta_b, spindle_b);
        % Save
        results_row = addNestedFieldsToResultsRow(complexity_measures, 'comp', results_row);



        %% Calculations only on clean epochs
        % group the clean epochs by sleep stage for further sleep stage specific analysis
        stage_epochs = group_epochs_by_stage(EEG_clean.clean_epochs);



        % Power spectrum analysis
        % Initialize structs to hold power spectrum results
        disp('- Extract power spectrum features');
        power_spectrum_results = struct('NREM1', [], 'NREM2', [], 'NREM3', [], 'REM', [], 'Overall', []);
        bands = {'delta', 'theta', 'alpha', 'beta', 'spindle'};
        regions = {'FC', 'PO'};  % Define two regions: Fronto-central (FC) and Parieto-occipital (PO)
        overall_absolute_power = struct('FC', [], 'PO', []);
        overall_relative_power = struct('FC', [], 'PO', []);
        stage_absolute_power = struct('FC', [], 'PO', []);
        stage_relative_power = struct('FC', [], 'PO', []);
        % Define channel groupings for each region
        region_channels.FC = find(contains({EEG_clean.chanlocs.labels}, {'F', 'C'}));  % Channels with labels F or C
        region_channels.PO = find(contains({EEG_clean.chanlocs.labels}, {'P', 'O'}));  % Channels with labels P or O
        % Initialize arrays to accumulate power values across all stages for the "Overall" calculation
        overall_absolute_power.FC = [];
        overall_relative_power.FC = [];
        overall_absolute_power.PO = [];
        overall_relative_power.PO = [];
        % Loop through each sleep stage and calculate power spectrum for each epoch
        sleep_stages = fieldnames(stage_epochs);  % NREM1, NREM2, NREM3, REM
        for stage_idx = 1:length(sleep_stages)
            stage_name = sleep_stages{stage_idx};  % Current stage (e.g., NREM1, NREM2, etc.)
            epochs = stage_epochs.(stage_name);    % Get epochs (start and end indices) for the current stage
            % Initialize arrays to hold power spectrum results for this stage and region
            stage_absolute_power.FC = [];
            stage_relative_power.FC = [];
            stage_absolute_power.PO = [];
            stage_relative_power.PO = [];
            % Loop over each region (FC and PO)
            for region_idx = 1:length(regions)
                region_name = regions{region_idx};
                region_chan = region_channels.(region_name);  % Get channel indices for this region
                % Loop over each epoch
                for epoch_idx = 1:size(epochs, 1)
                    start_idx = epochs(epoch_idx, 1);
                    end_idx = epochs(epoch_idx, 2);
                    indexs = start_idx:end_idx;  % Epoch indices
                    % Loop over each channel in the region
                    for chan = region_chan
                        % Calculate power spectrum for this channel and epoch
                        [absolute_power, relative_power, ~, ~] = calculate_power_spectrum(EEG_clean, chan, indexs, srate, freq_res, delta_b, theta_b, alpha_b, beta_b, spindle_b);
                        % Store the results for each epoch and channel in the stage-specific arrays
                        stage_absolute_power.(region_name) = [stage_absolute_power.(region_name); struct2array(absolute_power)];
                        stage_relative_power.(region_name) = [stage_relative_power.(region_name); struct2array(relative_power)];
                        % Accumulate results for the overall calculation
                        overall_absolute_power.(region_name) = [overall_absolute_power.(region_name); struct2array(absolute_power)];
                        overall_relative_power.(region_name) = [overall_relative_power.(region_name); struct2array(relative_power)];
                    end
                end
                % Store results in the struct for each sleep stage and region
                power_spectrum_results.(stage_name).(region_name).absolute_power = stage_absolute_power.(region_name);
                power_spectrum_results.(stage_name).(region_name).relative_power = stage_relative_power.(region_name);
                % Calculate mean and standard deviation for each band and region
                for band_idx = 1:length(bands)
                    band_name = bands{band_idx};
                    % Mean and SD of absolute power for this band in this stage and region
                    power_spectrum_results.(stage_name).(region_name).mean_absolute.(band_name) = mean(stage_absolute_power.(region_name)(:, band_idx), 'all');
                    power_spectrum_results.(stage_name).(region_name).std_absolute.(band_name) = std(stage_absolute_power.(region_name)(:, band_idx), 0, 'all');
                    % Mean and SD of relative power for this band in this stage and region
                    power_spectrum_results.(stage_name).(region_name).mean_relative.(band_name) = mean(stage_relative_power.(region_name)(:, band_idx), 'all');
                    power_spectrum_results.(stage_name).(region_name).std_relative.(band_name) = std(stage_relative_power.(region_name)(:, band_idx), 0, 'all');
                end
            end
        end
        % Now calculate the overall (whole night) values for each region
        for region_idx = 1:length(regions)
            region_name = regions{region_idx};
            % Overall mean and SD for each frequency band across all stages for this region
            for band_idx = 1:length(bands)
                band_name = bands{band_idx};
                % Mean and SD of absolute power for each band across all stages
                power_spectrum_results.Overall.(region_name).mean_absolute.(band_name) = mean(overall_absolute_power.(region_name)(:, band_idx), 'all');
                power_spectrum_results.Overall.(region_name).std_absolute.(band_name) = std(overall_absolute_power.(region_name)(:, band_idx), 0, 'all');
                % Mean and SD of relative power for each band across all stages
                power_spectrum_results.Overall.(region_name).mean_relative.(band_name) = mean(overall_relative_power.(region_name)(:, band_idx), 'all');
                power_spectrum_results.Overall.(region_name).std_relative.(band_name) = std(overall_relative_power.(region_name)(:, band_idx), 0, 'all');
            end
        end
        % Save
        results_row = addNestedFieldsToResultsRow(power_spectrum_results, 'power_spect', results_row);


        % FOOOF analysis
        % Channels grouped into FC and PO
        disp('- Extract fooof spectrum features');
        fc_channels = find(contains({EEG_clean.chanlocs.labels}, {'F', 'C'}));  % Find FC channels
        po_channels = find(contains({EEG_clean.chanlocs.labels}, {'P', 'O'}));  % Find PO channels
        % Initialize structs to hold FOOOF results for each stage and overall
        fooof_results = struct('NREM1', [], 'NREM2', [], 'NREM3', [], 'REM', [], 'Overall', []);
        bands = {'delta', 'theta', 'alpha', 'beta', 'spindle', 'gamma'};
        overall_normalized_absolute_power = struct('FC', [], 'PO', []);
        overall_peak_values = struct('FC', [], 'PO', []);
        stage_fooof_offset = struct('FC', [], 'PO', []);
        stage_fooof_exponent = struct('FC', [], 'PO', []);
        stage_normalized_absolute_power = struct('FC', [], 'PO', []);
        stage_peak_values = struct('FC', [], 'PO', []);
        % Initialize arrays to accumulate FOOOF values across all stages for the "Overall" calculation
        stage_fooof_offset.FC = [];
        stage_fooof_offset.PO = [];
        stage_fooof_exponent.FC = [];
        stage_fooof_exponent.PO = [];
        stage_normalized_absolute_power.FC = [];
        stage_normalized_absolute_power.PO = [];
        stage_peak_values.FC = [];
        stage_peak_values.PO = [];
        % Loop through each sleep stage and perform FOOOF analysis for each epoch
        sleep_stages = fieldnames(stage_epochs);  % NREM1, NREM2, NREM3, REM
        for stage_idx = 1:length(sleep_stages)
            stage_name = sleep_stages{stage_idx};  % Current stage (e.g., NREM1, NREM2, etc.)
            epochs = stage_epochs.(stage_name);    % Get epochs (start and end indices) for the current stage
            % Initialize arrays to hold FOOOF results for FC and PO channels separately
            stage_fooof_offset.FC = [];
            stage_fooof_offset.PO = [];
            stage_fooof_exponent.FC = [];
            stage_fooof_exponent.PO = [];
            stage_normalized_absolute_power.FC = [];
            stage_normalized_absolute_power.PO = [];
            stage_peak_values.FC = [];
            stage_peak_values.PO = [];
            for epoch_idx = 1:size(epochs, 1)
                start_idx = epochs(epoch_idx, 1);
                end_idx = epochs(epoch_idx, 2);
                indexs = start_idx:end_idx;  % Epoch indices
                % Process FC channels
                for chan = fc_channels
                    % Calculate power spectrum for this channel and epoch (needed for FOOOF)
                    [~, ~, freq, psd] = calculate_power_spectrum(EEG_clean, chan, indexs, srate, freq_res, delta_b, theta_b, alpha_b, beta_b, spindle_b);
                    % Run FOOOF and detect strongest peaks
                    [fooof_offset, fooof_exponent, peak_table, normalized_absolute_power] = run_fooof_and_detect_strongest_peaks(freq, psd, delta_b, theta_b, alpha_b, beta_b, spindle_b, gamma_b);
                    % Store the results for FC channels
                    stage_fooof_offset.FC = [stage_fooof_offset.FC; fooof_offset];
                    stage_fooof_exponent.FC = [stage_fooof_exponent.FC; fooof_exponent];
                    stage_normalized_absolute_power.FC = [stage_normalized_absolute_power.FC; struct2array(normalized_absolute_power)];
                    stage_peak_values.FC = [stage_peak_values.FC; struct2array(peak_table)];
                end
                % Process PO channels
                for chan = po_channels
                    % Calculate power spectrum for this channel and epoch (needed for FOOOF)
                    [~, ~, freq, psd] = calculate_power_spectrum(EEG_clean, chan, indexs, srate, freq_res, delta_b, theta_b, alpha_b, beta_b, spindle_b);
                    % Run FOOOF and detect strongest peaks
                    [fooof_offset, fooof_exponent, peak_table, normalized_absolute_power] = run_fooof_and_detect_strongest_peaks(freq, psd, delta_b, theta_b, alpha_b, beta_b, spindle_b, gamma_b);
                    % Store the results for PO channels
                    stage_fooof_offset.PO = [stage_fooof_offset.PO; fooof_offset];
                    stage_fooof_exponent.PO = [stage_fooof_exponent.PO; fooof_exponent];
                    stage_normalized_absolute_power.PO = [stage_normalized_absolute_power.PO; struct2array(normalized_absolute_power)];
                    stage_peak_values.PO = [stage_peak_values.PO; struct2array(peak_table)];
                end
            end
            % Store results in the struct for each sleep stage
            % For FC channels
            fooof_results.(stage_name).FC.fooof_offset = stage_fooof_offset.FC;
            fooof_results.(stage_name).FC.fooof_exponent = stage_fooof_exponent.FC;
            fooof_results.(stage_name).FC.normalized_absolute_power = stage_normalized_absolute_power.FC;
            fooof_results.(stage_name).FC.peak_values = stage_peak_values.FC;
            % For PO channels
            fooof_results.(stage_name).PO.fooof_offset = stage_fooof_offset.PO;
            fooof_results.(stage_name).PO.fooof_exponent = stage_fooof_exponent.PO;
            fooof_results.(stage_name).PO.normalized_absolute_power = stage_normalized_absolute_power.PO;
            fooof_results.(stage_name).PO.peak_values = stage_peak_values.PO;
            % Calculate mean and standard deviation for each parameter and overall
            % FC Channels
            fooof_results.(stage_name).FC.mean_fooof_offset = mean(stage_fooof_offset.FC, 'omitnan');
            fooof_results.(stage_name).FC.std_fooof_offset = std(stage_fooof_offset.FC, 0, 'omitnan');
            fooof_results.(stage_name).FC.mean_fooof_exponent = mean(stage_fooof_exponent.FC, 'omitnan');
            fooof_results.(stage_name).FC.std_fooof_exponent = std(stage_fooof_exponent.FC, 0, 'omitnan');
            % Mean and SD for normalized absolute power for each band for FC
            for band_idx = 1:length(bands)
                band_name = bands{band_idx};
                fooof_results.(stage_name).FC.mean_normalized_absolute_power.(band_name) = mean(stage_normalized_absolute_power.FC(:, band_idx), 'omitnan');
                fooof_results.(stage_name).FC.std_normalized_absolute_power.(band_name) = std(stage_normalized_absolute_power.FC(:, band_idx), 0, 'omitnan');
            end
            % Mean and SD for peak values (frequency, amplitude, bandwidth) for FC
            for field = fieldnames(peak_table)'  % Loop over all field names in peak_table
                field_name = field{1};  % Get the current field name as a string
                % Extract the column corresponding to the current field from the peak values
                peak_values_column = stage_peak_values.FC(:, ismember(fieldnames(peak_table), field_name));
                % Calculate the mean and standard deviation, ignoring NaNs
                fooof_results.(stage_name).FC.mean_peak_values.(field_name) = mean(peak_values_column, 'omitnan');
                fooof_results.(stage_name).FC.std_peak_values.(field_name) = std(peak_values_column, 0, 'omitnan');
            end
            % PO Channels (same process as above)
            fooof_results.(stage_name).PO.mean_fooof_offset = mean(stage_fooof_offset.PO, 'omitnan');
            fooof_results.(stage_name).PO.std_fooof_offset = std(stage_fooof_offset.PO, 0, 'omitnan');
            fooof_results.(stage_name).PO.mean_fooof_exponent = mean(stage_fooof_exponent.PO, 'omitnan');
            fooof_results.(stage_name).PO.std_fooof_exponent = std(stage_fooof_exponent.PO, 0, 'omitnan');
            for band_idx = 1:length(bands)
                band_name = bands{band_idx};
                fooof_results.(stage_name).PO.mean_normalized_absolute_power.(band_name) = mean(stage_normalized_absolute_power.PO(:, band_idx), 'omitnan');
                fooof_results.(stage_name).PO.std_normalized_absolute_power.(band_name) = std(stage_normalized_absolute_power.PO(:, band_idx), 0, 'omitnan');
            end
            for field = fieldnames(peak_table)'
                field_name = field{1};
                peak_values_column = stage_peak_values.PO(:, ismember(fieldnames(peak_table), field_name));
                fooof_results.(stage_name).PO.mean_peak_values.(field_name) = mean(peak_values_column, 'omitnan');
                fooof_results.(stage_name).PO.std_peak_values.(field_name) = std(peak_values_column, 0, 'omitnan');
            end
        end
        % Save
        results_row = addNestedFieldsToResultsRow(fooof_results, 'fooof', results_row);



        % Run the slow wave detection in YASA
        disp('- Extract SW features');
        coupling_py = py.dict(pyargs('freq_sp', 12:16, 'p', 0.05, 'time', 1));
        sws_detection_results = py.yasa.sw_detect(data = clean_data_for_py, ...
            sf = py.float(EEG_clean.srate), ...
            ch_names = py.list({EEG_clean.chanlocs.labels}), ...
            hypno = py.numpy.array(hypnogram_clean_upsampled, 'int'), ...
            include = py.numpy.array([2, 3], 'int'), ...  % stages to include in SWS detection (NREM2 and NREM3)
            freq_sw = py.tuple([0.5, 4.5]), ...
            dur_neg = py.tuple([0.1, 1.5]), ...
            dur_pos = py.tuple([0.1, 1.5]), ...
            amp_neg = py.tuple([5, 1000]), ...
            amp_pos = py.tuple([5, 1000]), ...
            amp_ptp = py.tuple([20, 1000]), ...
            coupling = coupling_py, ...  % Disable phase coupling calculation
            remove_outliers = false, ...  % Disable outlier removal
            verbose = false);
        % Make the python structure accessible for MATLAB
        if ~isempty(sws_detection_results)
            sws_detection_info = sws_detection_results.summary();
            sws_colnames = cellfun(@string, cell(sws_detection_info.columns.values.tolist()));
            sws_detection_info = cell(sws_detection_info.values.tolist());
            sws_detection_info = cellfun(@(xs)cell2table(cell(xs), 'VariableNames', sws_colnames), sws_detection_info, 'UniformOutput', false);
            sws_detection_info = vertcat(sws_detection_info{:});
            sws_detection_info.Stage = cellfun(@double, sws_detection_info.Stage);
            sws_detection_info.Channel = string(sws_detection_info.Channel);
            sws_detection_info.IdxChannel = cellfun(@double, sws_detection_info.IdxChannel);
            sws_detection_info.neg_sw_peaks_ind = sws_detection_info.NegPeak * EEG_clean.srate;
            sws_detection_info.duration_negative_phase = (sws_detection_info.MidCrossing - sws_detection_info.Start) * EEG_clean.srate;
            sws_detection_info.duration_positive_phase = (sws_detection_info.End - sws_detection_info.MidCrossing) * EEG_clean.srate;
            sws_detection_info.duration_start_to_neg = (sws_detection_info.NegPeak - sws_detection_info.Start) * EEG_clean.srate;
            sws_detection_info.duration_neg_to_mid = (sws_detection_info.MidCrossing - sws_detection_info.NegPeak) * EEG_clean.srate;
            sws_detection_info.duration_mid_to_pos = (sws_detection_info.PosPeak - sws_detection_info.MidCrossing) * EEG_clean.srate;
            sws_detection_info.duration_pos_to_end = (sws_detection_info.End - sws_detection_info.PosPeak) * EEG_clean.srate;
        end
        % Initialize an empty struct to hold the results for each cycle
        SW_results = struct();
        % Loop through each sleep cycle in the cycle_table
        for cycle_idx = 1:height(cycle_table)
            start_epoch = cycle_table.start_epoch(cycle_idx) * 30;  % Convert to data points
            end_epoch = cycle_table.end_epoch(cycle_idx) * 30;      % Convert to data points
            % Initialize arrays to hold slow wave data for this cycle (NREM stages only)
            sw_cycle_data_NREM = [];
            % Loop through the detected slow waves in sws_detection_info
            for sw_idx = 1:height(sws_detection_info)
                sw_start = sws_detection_info.Start(sw_idx);  % Start of slow wave (in data points)
                sw_stage = sws_detection_info.Stage(sw_idx);  % Stage of the slow wave
                % Check if the slow wave occurs within the current cycle and is in NREM (stage 2 or 3)
                if sw_start >= start_epoch && sw_start <= end_epoch && (sw_stage == 2 || sw_stage == 3)
                    sw_cycle_data_NREM = [sw_cycle_data_NREM; sws_detection_info(sw_idx, :)];
                end
            end
            % For NREM slow waves, calculate mean and std of the relevant metrics
            if ~isempty(sw_cycle_data_NREM)
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration = mean(sw_cycle_data_NREM.Duration);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration = std(sw_cycle_data_NREM.Duration);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ValNegPeak = mean(sw_cycle_data_NREM.ValNegPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ValNegPeak = std(sw_cycle_data_NREM.ValNegPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ValPosPeak = mean(sw_cycle_data_NREM.ValPosPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ValPosPeak = std(sw_cycle_data_NREM.ValPosPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.PTP = mean(sw_cycle_data_NREM.PTP);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.PTP = std(sw_cycle_data_NREM.PTP);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.Slope = mean(sw_cycle_data_NREM.Slope);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.Slope = std(sw_cycle_data_NREM.Slope);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.Frequency = mean(sw_cycle_data_NREM.Frequency);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.Frequency = std(sw_cycle_data_NREM.Frequency);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.SigmaPeak = mean(sw_cycle_data_NREM.SigmaPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.SigmaPeak = std(sw_cycle_data_NREM.SigmaPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.PhaseAtSigmaPeak = mean(sw_cycle_data_NREM.PhaseAtSigmaPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.PhaseAtSigmaPeak = std(sw_cycle_data_NREM.PhaseAtSigmaPeak);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ndPAC = mean(sw_cycle_data_NREM.ndPAC);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ndPAC = std(sw_cycle_data_NREM.ndPAC);
                % Durations
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration_negative_phase = mean(sw_cycle_data_NREM.duration_negative_phase);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration_negative_phase = std(sw_cycle_data_NREM.duration_negative_phase);
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration_positive_phase = mean(sw_cycle_data_NREM.duration_positive_phase);
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration_positive_phase = std(sw_cycle_data_NREM.duration_positive_phase);
            else
                % If no NREM slow waves found in this cycle, set fields to NaN
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ValNegPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ValNegPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ValPosPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ValPosPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.PTP = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.PTP = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.Slope = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.Slope = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.Frequency = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.Frequency = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.SigmaPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.SigmaPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.PhaseAtSigmaPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.PhaseAtSigmaPeak = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.ndPAC = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.ndPAC = NaN;
                % Durations
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration_negative_phase = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration_negative_phase = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).mean.duration_positive_phase = NaN;
                SW_results.(['Cycle_' num2str(cycle_idx)]).std.duration_positive_phase = NaN;
            end
        end
        % Calculate the difference between the first and last cycle
        first_cycle = SW_results.Cycle_1;
        % Get the last cycle data
        % Loop through all cycles from the end to find the last valid one
        for cycle_idx = height(cycle_table):-1:1
            % Access the current cycle
            current_cycle = SW_results.(['Cycle_' num2str(cycle_idx)]);
            % Check if mean is a struct and has fields
            if isstruct(current_cycle.mean) && ~isempty(fieldnames(current_cycle.mean))
                % Check if all entries in current_cycle.mean are NaN
                if ~all(structfun(@(x) all(isnan(x)), current_cycle.mean))
                    % If any entry is not NaN, assign it to last_cycle and break the loop
                    last_cycle = current_cycle;
                    break;
                end
            else
                % If mean is NaN, skip to the next cycle
                continue;
            end
        end
        % Calculate the difference in mean values between the first and last cycle
        SW_results.Difference_first_last_cycle.mean.duration = first_cycle.mean.duration - last_cycle.mean.duration;
        SW_results.Difference_first_last_cycle.mean.ValNegPeak = first_cycle.mean.ValNegPeak - last_cycle.mean.ValNegPeak;
        SW_results.Difference_first_last_cycle.mean.ValPosPeak = first_cycle.mean.ValPosPeak - last_cycle.mean.ValPosPeak;
        SW_results.Difference_first_last_cycle.mean.PTP = first_cycle.mean.PTP - last_cycle.mean.PTP;
        SW_results.Difference_first_last_cycle.mean.Slope = first_cycle.mean.Slope - last_cycle.mean.Slope;
        SW_results.Difference_first_last_cycle.mean.Frequency = first_cycle.mean.Frequency - last_cycle.mean.Frequency;
        SW_results.Difference_first_last_cycle.mean.SigmaPeak = first_cycle.mean.SigmaPeak - last_cycle.mean.SigmaPeak;
        SW_results.Difference_first_last_cycle.mean.PhaseAtSigmaPeak = first_cycle.mean.PhaseAtSigmaPeak - last_cycle.mean.PhaseAtSigmaPeak;
        SW_results.Difference_first_last_cycle.mean.ndPAC = first_cycle.mean.ndPAC - last_cycle.mean.ndPAC;
        SW_results.Difference_first_last_cycle.mean.duration_negative_phase = first_cycle.mean.duration_negative_phase - last_cycle.mean.duration_negative_phase;
        SW_results.Difference_first_last_cycle.mean.duration_positive_phase = first_cycle.mean.duration_positive_phase - last_cycle.mean.duration_positive_phase;
        % Calculate the difference in standard deviation between the first and last cycle
        SW_results.Difference_first_last_cycle.std.duration = first_cycle.std.duration - last_cycle.std.duration;
        SW_results.Difference_first_last_cycle.std.ValNegPeak = first_cycle.std.ValNegPeak - last_cycle.std.ValNegPeak;
        SW_results.Difference_first_last_cycle.std.ValPosPeak = first_cycle.std.ValPosPeak - last_cycle.std.ValPosPeak;
        SW_results.Difference_first_last_cycle.std.PTP = first_cycle.std.PTP - last_cycle.std.PTP;
        SW_results.Difference_first_last_cycle.std.Slope = first_cycle.std.Slope - last_cycle.std.Slope;
        SW_results.Difference_first_last_cycle.std.Frequency = first_cycle.std.Frequency - last_cycle.std.Frequency;
        SW_results.Difference_first_last_cycle.std.SigmaPeak = first_cycle.std.SigmaPeak - last_cycle.std.SigmaPeak;
        SW_results.Difference_first_last_cycle.std.PhaseAtSigmaPeak = first_cycle.std.PhaseAtSigmaPeak - last_cycle.std.PhaseAtSigmaPeak;
        SW_results.Difference_first_last_cycle.std.ndPAC = first_cycle.std.ndPAC - last_cycle.std.ndPAC;
        SW_results.Difference_first_last_cycle.std.duration_negative_phase = first_cycle.std.duration_negative_phase - last_cycle.std.duration_negative_phase;
        SW_results.Difference_first_last_cycle.std.duration_positive_phase = first_cycle.std.duration_positive_phase - last_cycle.std.duration_positive_phase;
        % Initialize empty matrices to store all cycle data for overall calculation
        all_cycle_data_mean = [];
        all_cycle_data_std = [];
        % Loop through each cycle and aggregate the data
        for cycle_idx = 1:height(cycle_table)
            cycle_name = ['Cycle_' num2str(cycle_idx)];
            if isfield(SW_results, cycle_name)
                % Check if 'mean' is a struct with fields, indicating a valid cycle
                if isstruct(SW_results.(cycle_name).mean)
                    % Extract and concatenate the data for overall calculation, ignoring NaNs
                    mean_values = struct2array(SW_results.(cycle_name).mean);
                    std_values = struct2array(SW_results.(cycle_name).std);
                    % Append non-NaN data to the overall data matrices
                    if ~isempty(mean_values)
                        all_cycle_data_mean = [all_cycle_data_mean; mean_values];
                    end
                    if ~isempty(std_values)
                        all_cycle_data_std = [all_cycle_data_std; std_values];
                    end
                end
            end
        end
        % Calculate the overall mean across all cycles
        if ~isempty(all_cycle_data_mean)
            overall_mean_values = mean(all_cycle_data_mean, 1, 'omitnan');
            overall_std_values = std(all_cycle_data_mean, 0, 1, 'omitnan');
            % Store in the SW_results structure
            SW_results.Overall.mean.duration = overall_mean_values(1);
            SW_results.Overall.mean.ValNegPeak = overall_mean_values(2);
            SW_results.Overall.mean.ValPosPeak = overall_mean_values(3);
            SW_results.Overall.mean.PTP = overall_mean_values(4);
            SW_results.Overall.mean.Slope = overall_mean_values(5);
            SW_results.Overall.mean.Frequency = overall_mean_values(6);
            SW_results.Overall.mean.SigmaPeak = overall_mean_values(7);
            SW_results.Overall.mean.PhaseAtSigmaPeak = overall_mean_values(8);
            SW_results.Overall.mean.ndPAC = overall_mean_values(9);
            SW_results.Overall.mean.duration_negative_phase = overall_mean_values(10);
            SW_results.Overall.mean.duration_positive_phase = overall_mean_values(11);
            % Store the overall standard deviation
            SW_results.Overall.std.duration = overall_std_values(1);
            SW_results.Overall.std.ValNegPeak = overall_std_values(2);
            SW_results.Overall.std.ValPosPeak = overall_std_values(3);
            SW_results.Overall.std.PTP = overall_std_values(4);
            SW_results.Overall.std.Slope = overall_std_values(5);
            SW_results.Overall.std.Frequency = overall_std_values(6);
            SW_results.Overall.std.SigmaPeak = overall_std_values(7);
            SW_results.Overall.std.PhaseAtSigmaPeak = overall_std_values(8);
            SW_results.Overall.std.ndPAC = overall_std_values(9);
            SW_results.Overall.std.duration_negative_phase = overall_std_values(10);
            SW_results.Overall.std.duration_positive_phase = overall_std_values(11);
        end
        % Save
        results_row = addNestedFieldsToResultsRow(SW_results.Difference_first_last_cycle, 'SW_diff_first_last', results_row);
        results_row = addNestedFieldsToResultsRow(SW_results.Overall, 'SW_overall', results_row);



        % Run the spindle detection from YASA
        % Correct the 'thresh' parameter to be a Python dictionary
        disp('- Extract spindle features');
        thresh_py = py.dict(pyargs('corr', 0.65, 'rel_pow', 0.2, 'rms', 1.5));
        spindles_detection_results = py.yasa.spindles_detect(data = clean_data_for_py, ...
            sf = py.float(EEG_clean.srate), ...
            ch_names = py.list({EEG_clean.chanlocs.labels}), ...
            hypno = py.numpy.array(hypnogram_clean_upsampled, 'int'), ...
            include = py.numpy.array([2, 3], 'int'), ...
            freq_sp = py.tuple([12, 15]), ...
            freq_broad = py.tuple([1, 30]), ...
            duration = py.tuple([0.5, 2]), ...
            min_distance = py.int(500), ...
            thresh = thresh_py, ...
            multi_only = py.False, ...
            remove_outliers = py.False, ...
            verbose = py.False);
        % Extract the spindle summary from the detection results
        if ~isempty(spindles_detection_results)
            spindles_summary_py = spindles_detection_results.summary();
            % Extract the column names
            spindle_colnames = cellfun(@string, cell(spindles_summary_py.columns.values.tolist()));
            % Convert the Python dataframe into a MATLAB cell array
            spindles_detection_info = cell(spindles_summary_py.values.tolist());
            % Convert each row in the cell array to a table with the extracted column names
            spindles_detection_info = cellfun(@(xs) cell2table(cell(xs), 'VariableNames', spindle_colnames), spindles_detection_info, 'UniformOutput', false);
            % Concatenate all rows to form the final table
            spindles_detection_info = vertcat(spindles_detection_info{:});
            % Convert necessary columns to appropriate data types
            spindles_detection_info.Stage = cellfun(@double, spindles_detection_info.Stage);  % Convert 'Stage' to double
            spindles_detection_info.Channel = string(spindles_detection_info.Channel);        % Convert 'Channel' to string
            spindles_detection_info.IdxChannel = cellfun(@double, spindles_detection_info.IdxChannel);  % Convert 'IdxChannel' to double
        end
        % Initialize an empty struct to hold the spindle results for each cycle
        Spindle_results = struct();
        % Loop through each sleep cycle in the cycle_table
        for cycle_idx = 1:height(cycle_table)
            start_epoch = cycle_table.start_epoch(cycle_idx) * 30;  % Convert to data points
            end_epoch = cycle_table.end_epoch(cycle_idx) * 30;      % Convert to data points
            % Initialize arrays to hold spindle data for this cycle (NREM stages only)
            spindle_cycle_data_NREM = [];
            % Loop through the detected spindles in spindles_detection_info
            for spindle_idx = 1:height(spindles_detection_info)
                spindle_start = spindles_detection_info.Start(spindle_idx);  % Start of spindle (in data points)
                spindle_stage = spindles_detection_info.Stage(spindle_idx);  % Stage of the spindle
                % Check if the spindle occurs within the current cycle and is in NREM (stage 2 or 3)
                if spindle_start >= start_epoch && spindle_start <= end_epoch && (spindle_stage == 2 || spindle_stage == 3)
                    spindle_cycle_data_NREM = [spindle_cycle_data_NREM; spindles_detection_info(spindle_idx, :)];
                end
            end
            % For NREM spindles, calculate mean and std of the relevant metrics
            if ~isempty(spindle_cycle_data_NREM)
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.duration = mean(spindle_cycle_data_NREM.Duration);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.duration = std(spindle_cycle_data_NREM.Duration);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Amplitude = mean(spindle_cycle_data_NREM.Amplitude);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Amplitude = std(spindle_cycle_data_NREM.Amplitude);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.RMS = mean(spindle_cycle_data_NREM.RMS);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.RMS = std(spindle_cycle_data_NREM.RMS);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.AbsPower = mean(spindle_cycle_data_NREM.AbsPower);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.AbsPower = std(spindle_cycle_data_NREM.AbsPower);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.RelPower = mean(spindle_cycle_data_NREM.RelPower);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.RelPower = std(spindle_cycle_data_NREM.RelPower);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Frequency = mean(spindle_cycle_data_NREM.Frequency);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Frequency = std(spindle_cycle_data_NREM.Frequency);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Oscillations = mean(spindle_cycle_data_NREM.Oscillations);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Oscillations = std(spindle_cycle_data_NREM.Oscillations);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Symmetry = mean(spindle_cycle_data_NREM.Symmetry);
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Symmetry = std(spindle_cycle_data_NREM.Symmetry);
            else
                % If no NREM spindles found in this cycle, set fields to NaN
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.duration = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.duration = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Amplitude = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Amplitude = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.RMS = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.RMS = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.AbsPower = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.AbsPower = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.RelPower = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.RelPower = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Frequency = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Frequency = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Oscillations = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Oscillations = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean.Symmetry = NaN;
                Spindle_results.(['Cycle_' num2str(cycle_idx)]).std.Symmetry = NaN;
            end
        end
        % Calculate the difference between the first and last cycle
        first_cycle = Spindle_results.Cycle_1;
        % Loop through all cycles from the end to find the last valid one
        for cycle_idx = height(cycle_table):-1:1
            % Access the current cycle
            current_cycle = Spindle_results.(['Cycle_' num2str(cycle_idx)]);
            % Check if mean is a struct and has fields
            if isstruct(current_cycle.mean) && ~isempty(fieldnames(current_cycle.mean))
                % Check if all entries in current_cycle.mean are NaN
                if ~all(structfun(@(x) all(isnan(x)), current_cycle.mean))
                    % If any entry is not NaN, assign it to last_cycle and break the loop
                    last_cycle = current_cycle;
                    break;
                end
            else
                % If mean is NaN, skip to the next cycle
                continue;
            end
        end
        % Calculate the difference in mean values between the first and last cycle
        Spindle_results.Difference_first_last_cycle.mean.duration = first_cycle.mean.duration - last_cycle.mean.duration;
        Spindle_results.Difference_first_last_cycle.mean.Amplitude = first_cycle.mean.Amplitude - last_cycle.mean.Amplitude;
        Spindle_results.Difference_first_last_cycle.mean.RMS = first_cycle.mean.RMS - last_cycle.mean.RMS;
        Spindle_results.Difference_first_last_cycle.mean.AbsPower = first_cycle.mean.AbsPower - last_cycle.mean.AbsPower;
        Spindle_results.Difference_first_last_cycle.mean.RelPower = first_cycle.mean.RelPower - last_cycle.mean.RelPower;
        Spindle_results.Difference_first_last_cycle.mean.Frequency = first_cycle.mean.Frequency - last_cycle.mean.Frequency;
        Spindle_results.Difference_first_last_cycle.mean.Oscillations = first_cycle.mean.Oscillations - last_cycle.mean.Oscillations;
        Spindle_results.Difference_first_last_cycle.mean.Symmetry = first_cycle.mean.Symmetry - last_cycle.mean.Symmetry;
        % Calculate the difference in standard deviation between the first and last cycle
        Spindle_results.Difference_first_last_cycle.std.duration = first_cycle.std.duration - last_cycle.std.duration;
        Spindle_results.Difference_first_last_cycle.std.Amplitude = first_cycle.std.Amplitude - last_cycle.std.Amplitude;
        Spindle_results.Difference_first_last_cycle.std.RMS = first_cycle.std.RMS - last_cycle.std.RMS;
        Spindle_results.Difference_first_last_cycle.std.AbsPower = first_cycle.std.AbsPower - last_cycle.std.AbsPower;
        Spindle_results.Difference_first_last_cycle.std.RelPower = first_cycle.std.RelPower - last_cycle.std.RelPower;
        Spindle_results.Difference_first_last_cycle.std.Frequency = first_cycle.std.Frequency - last_cycle.std.Frequency;
        Spindle_results.Difference_first_last_cycle.std.Oscillations = first_cycle.std.Oscillations - last_cycle.std.Oscillations;
        Spindle_results.Difference_first_last_cycle.std.Symmetry = first_cycle.std.Symmetry - last_cycle.std.Symmetry;
        % Initialize empty matrices to store all cycle data for overall calculation
        all_cycle_data_mean = [];
        all_cycle_data_std = [];
         for cycle_idx = 1:height(cycle_table)
            cycle_name = ['Cycle_' num2str(cycle_idx)];
            if isfield(Spindle_results, cycle_name)
                % Check if 'mean' is a struct with fields, indicating a valid cycle
                if isstruct(Spindle_results.(cycle_name).mean)
                    % Extract and concatenate the data for overall calculation, ignoring NaNs
                    mean_values = struct2array(Spindle_results.(cycle_name).mean);
                    std_values = struct2array(Spindle_results.(cycle_name).std);
                    % Append non-NaN data to the overall data matrices
                    if ~isempty(mean_values)
                        all_cycle_data_mean = [all_cycle_data_mean; mean_values];
                    end
                    if ~isempty(std_values)
                        all_cycle_data_std = [all_cycle_data_std; std_values];
                    end
                end
            end
        end
        % Calculate the overall mean and standard deviation across all cycles
        if ~isempty(all_cycle_data_mean)
            overall_mean_values = mean(all_cycle_data_mean, 1, 'omitnan');
            overall_std_values = std(all_cycle_data_mean, 0, 1, 'omitnan');
            % Store in the Spindle_results structure
            Spindle_results.Overall.mean.duration = overall_mean_values(1);
            Spindle_results.Overall.mean.Amplitude = overall_mean_values(2);
            Spindle_results.Overall.mean.RMS = overall_mean_values(3);
            Spindle_results.Overall.mean.AbsPower = overall_mean_values(4);
            Spindle_results.Overall.mean.RelPower = overall_mean_values(5);
            Spindle_results.Overall.mean.Frequency = overall_mean_values(6);
            Spindle_results.Overall.mean.Oscillations = overall_mean_values(7);
            Spindle_results.Overall.mean.Symmetry = overall_mean_values(8);
            % Store the overall standard deviation
            Spindle_results.Overall.std.duration = overall_std_values(1);
            Spindle_results.Overall.std.Amplitude = overall_std_values(2);
            Spindle_results.Overall.std.RMS = overall_std_values(3);
            Spindle_results.Overall.std.AbsPower = overall_std_values(4);
            Spindle_results.Overall.std.RelPower = overall_std_values(5);
            Spindle_results.Overall.std.Frequency = overall_std_values(6);
            Spindle_results.Overall.std.Oscillations = overall_std_values(7);
            Spindle_results.Overall.std.Symmetry = overall_std_values(8);
        end
        % Save
        results_row = addNestedFieldsToResultsRow(Spindle_results.Difference_first_last_cycle, 'Spindle_diff_first_last', results_row);
        results_row = addNestedFieldsToResultsRow(Spindle_results.Overall, 'Spindle_overall', results_row);



        % Detect phasic and tonic REM
        % Define parameters
        disp('- Detect phasic and tonic REM');
        fs = EEG_clean.srate;  % Sampling rate
        eog_signal = EEG_clean.data_EOG(1, :);  % Single EOG channel
        freq_range = [0.5, 5];  % Frequency range for REM eye movements
        epoch_duration = 30;  % Duration of each epoch in seconds
        % Initialize storage for results
        rem_analysis = struct();
        total_phasic_duration = 0;
        total_tonic_duration = 0;
        % Loop over each sleep cycle
        for cycle_idx = 1:height(cycle_table)
            % Define the start and end points for the current cycle in terms of samples
            start_epoch = cycle_table.start_epoch(cycle_idx);
            end_epoch = cycle_table.end_epoch(cycle_idx);
            start_sample = (start_epoch - 1) * epoch_duration * fs + 1;
            end_sample = end_epoch * epoch_duration * fs;
            % Extract the EOG signal for the current cycle
            eog_cycle_signal = eog_signal(start_sample:end_sample);
            % Mask to select only REM segments within the current cycle
            cycle_rem_mask = hypnogram_upsampled(start_sample:end_sample) == 5;
            eog_signal_rem = eog_cycle_signal(cycle_rem_mask);  % Use only REM data for this cycle
            % Bandpass filter the REM EOG signal for the cycle
            eog_filtered_rem = bandpass(eog_signal_rem, freq_range, fs);
            % Compute the amplitude envelope of the REM EOG signal
            eog_envelope_rem = abs(hilbert(eog_filtered_rem));
            % Calculate an amplitude threshold only on the REM data
            amp_threshold = 1.5 * mean(eog_envelope_rem);
            % Detect phasic REM within REM segments where amplitude is above the threshold
            phasic_mask_rem = eog_envelope_rem > amp_threshold;
            % Calculate phasic and tonic REM durations based on REM-only data
            phasic_duration = sum(phasic_mask_rem) / fs;
            tonic_duration = (length(phasic_mask_rem) - sum(phasic_mask_rem)) / fs;
            % Total REM duration for this cycle
            total_rem_duration = phasic_duration + tonic_duration;
            % Update overall durations for phasic and tonic REM
            total_phasic_duration = total_phasic_duration + phasic_duration;
            total_tonic_duration = total_tonic_duration + tonic_duration;
            % Calculate phasic-to-tonic ratio for this cycle
            phasic_to_tonic_ratio = phasic_duration / tonic_duration;
            % Store results for the current cycle in a structured format
            rem_analysis.(['Cycle_' num2str(cycle_idx)]).phasic_duration_min = phasic_duration / 60;
            rem_analysis.(['Cycle_' num2str(cycle_idx)]).tonic_duration_min = tonic_duration / 60;
            rem_analysis.(['Cycle_' num2str(cycle_idx)]).total_rem_duration_min = total_rem_duration / 60;
            rem_analysis.(['Cycle_' num2str(cycle_idx)]).phasic_to_tonic_ratio = phasic_to_tonic_ratio;
        end
        % Calculate overall phasic and tonic durations and ratio
        overall_phasic_to_tonic_ratio = total_phasic_duration / total_tonic_duration;
        % Store overall results in the structure
        rem_analysis.Overall.total_phasic_duration_min = total_phasic_duration / 60;
        rem_analysis.Overall.total_tonic_duration_min = total_tonic_duration / 60;
        rem_analysis.Overall.phasic_to_tonic_ratio = overall_phasic_to_tonic_ratio;
        % Calculate difference values based on the first and last cycle
        if height(cycle_table) > 1
            % Initialize variables to store the last valid values
            last_cycle_ratio = NaN;
            last_cycle_rem_duration = NaN;
            % Loop from the end to find the last valid cycle where both values are not NaN
            for cycle_idx = height(cycle_table):-1:1
                % Access the current cycle
                current_cycle = rem_analysis.(['Cycle_' num2str(cycle_idx)]);
                % Check if both fields are not NaN
                if ~isnan(current_cycle.phasic_to_tonic_ratio) && ~isnan(current_cycle.total_rem_duration_min)
                    last_cycle_ratio = current_cycle.phasic_to_tonic_ratio;
                    last_cycle_rem_duration = current_cycle.total_rem_duration_min;
                    break;
                end
            end
            % Calculate the phasic-to-tonic ratio difference with the first cycle
            first_cycle_ratio = rem_analysis.Cycle_1.phasic_to_tonic_ratio;
            rem_analysis.Difference_first_last_cycle.phasic_to_tonic_ratio_ratio = last_cycle_ratio / first_cycle_ratio;
            % Similarly, calculate the total REM duration ratio
            first_cycle_rem_duration = rem_analysis.Cycle_1.total_rem_duration_min;
            rem_analysis.Difference_first_last_cycle.total_rem_duration_ratio = last_cycle_rem_duration / first_cycle_rem_duration;
        else
            warning('Only one cycle detected; cannot calculate difference between first and last cycle.');
            rem_analysis.Difference_first_last_cycle.phasic_to_tonic_ratio_diff = NaN;
            rem_analysis.Difference_first_last_cycle.total_rem_duration_ratio = NaN;
        end
        % Save
        results_row = addNestedFieldsToResultsRow(rem_analysis.Difference_first_last_cycle, 'REM_diff_first_last', results_row);
        results_row = addNestedFieldsToResultsRow(rem_analysis.Overall, 'REM_overall', results_row);




        % % Plot amplitude distribution
        % % Define the amplitude threshold for phasic REM
        % amp_threshold = 1.5 * mean(eog_envelope_rem);  % Example threshold calculation
        % % Plot the histogram of the amplitude envelope
        % figure;
        % histogram(eog_envelope_rem, 'BinWidth', 0.5);  % Adjust BinWidth as needed
        % hold on;
        % % Plot the cutoff line
        % y_limits = ylim;  % Get the current y-axis limits
        % plot([amp_threshold, amp_threshold], y_limits, 'r--', 'LineWidth', 2);
        % % Add title and labels
        % title('Histogram of Amplitude Envelope in REM EOG Signal');
        % xlabel('Amplitude');
        % ylabel('Count');
        % % Add legend
        % legend('Amplitude', 'Cutoff for Phasic REM');
        % % Release hold
        % hold off;
        %
        %
        % % Plot timedistribution tonic phasic
        % % Define the time vector for the REM-only signal
        % time_rem = (1:length(eog_envelope_rem)) / fs;
        % % Define the amplitude threshold
        % amp_threshold = 1.5 * mean(eog_envelope_rem);  % Example threshold calculation
        % % Create binary data (1 for phasic, 0 for tonic)
        % phasic_tonic_binary = eog_envelope_rem > amp_threshold;
        % % Define the indices for the 4 to 6-minute interval
        % start_idx = round(240 * fs);  % 4 minutes in samples
        % end_idx = min(round(300 * fs), length(phasic_tonic_binary));  % 6 minutes in samples, limited by data length
        % % Extract the data for this interval
        % time_interval = time_rem(start_idx:end_idx);
        % phasic_tonic_interval = phasic_tonic_binary(start_idx:end_idx);
        % % Plot the binary hypnogram-like data for the 4-6 min interval
        % figure;
        % plot(time_interval, phasic_tonic_interval, 'LineWidth', 1.5);  % Line plot for the selected interval
        % % Set axis limits and labels
        % ylim([-0.1, 1.1]);  % Limits for binary plot
        % yticks([0 1]);      % Set y-axis ticks to 0 and 1
        % yticklabels({'Tonic REM', 'Phasic REM'});
        % xlabel('Time (s)');
        % title('1 min Window of Phasic and Tonic REM Segments');
        % % Add a grid for clarity
        % grid on;



        % Calculate Phase Locking Value, Modulation Index and Mean Vektor Length for SW-Spindle-Coupling
        % PLV is a measure that quantifies the consistency of the phase relationship between two signals over time
        % MI tells you about the strength and specificity of the coupling
        % MVL provides insight into the directionality or phase preference of the amplitude modulation
        % Initialize the main structure to store PAC results for each channel
        disp('- Calculate Phase Locking values');
        pac_results = struct();
        n_channels = size(EEG_clean.clean_data, 1);  % Number of channels in the cleaned data
        % Initialize arrays to store overall values for averaging
        plv_values = zeros(n_channels, 1);
        mi_values = zeros(n_channels, 1);
        mvl_values = zeros(n_channels, 1);
        for chan = 1:n_channels
            % Get the channel name
            chan_name = EEG_clean.chanlocs(chan).labels;
            % Filter for low and high frequency bands
            low_freq_signal = bandpass(EEG_clean.clean_data(chan, :), delta_b, EEG_clean.srate);
            high_freq_signal = bandpass(EEG_clean.clean_data(chan, :), spindle_b, EEG_clean.srate);
            % Compute the phase of the low-frequency signal
            phase_low = angle(hilbert(low_freq_signal));
            % Compute the amplitude envelope for the high-frequency signal
            amp_env_high = abs(hilbert(high_freq_signal));
            phase_high = angle(hilbert(amp_env_high));
            % Calculate the Phase Locking Value (PLV)
            phase_diff = phase_low - phase_high;
            plv = abs(mean(exp(1i * phase_diff)));
            plv_values(chan) = plv;  % Store for overall calculation
            % Bin the phases and calculate the mean amplitude in each bin
            n_bins = 18;
            bin_edges = linspace(-pi, pi, n_bins + 1);
            amp_binned = zeros(1, n_bins);
            for bin = 1:n_bins
                idx = phase_low >= bin_edges(bin) & phase_low < bin_edges(bin + 1);
                amp_binned(bin) = mean(amp_env_high(idx));
            end
            % Normalize amplitude bins for Modulation Index (MI)
            p = amp_binned / sum(amp_binned);
            mi = (log(n_bins) + sum(p .* log(p))) / log(n_bins);
            mi_values(chan) = mi;  % Store for overall calculation
            % Calculate Mean Vector Length (MVL)
            phase_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;
            mvl = abs(sum(amp_binned .* exp(1i * phase_centers)) / sum(amp_binned));
            mvl_values(chan) = mvl;  % Store for overall calculation
            % Store essential results in the structure for the specific channel
            pac_results.(chan_name).plv = plv;
            pac_results.(chan_name).mi = mi;
            pac_results.(chan_name).mvl = mvl;
        end
        % Compute overall (mean) values across all channels
        pac_results.overall.plv = mean(plv_values);
        pac_results.overall.mi = mean(mi_values);
        pac_results.overall.mvl = mean(mvl_values);
        % Save
        results_row = addNestedFieldsToResultsRow(pac_results, 'PAC', results_row);


         
        % Convert results_row to table
        new_row_table = struct2table(results_row, 'AsArray', true);
        % Merge new_row_table into results_table
        if isempty(results_table)
            results_table = new_row_table;
        else
            % Combine variable names
            all_var_names = unique([results_table.Properties.VariableNames, new_row_table.Properties.VariableNames]);
            % Ensure both tables have all variables
            for varIdx = 1:length(all_var_names)
                varName = all_var_names{varIdx};
                if ~ismember(varName, results_table.Properties.VariableNames)
                    % Add missing variable to results_table
                    results_table.(varName) = repmat({NaN}, height(results_table), 1);
                end
                if ~ismember(varName, new_row_table.Properties.VariableNames)
                    % Add missing variable to new_row_table
                    new_row_table.(varName) = NaN;
                end
            end
            % Now append new_row_table to results_table
            results_table = [results_table; new_row_table];
        end
        % Save the updated results_table
        save(results_table_file, 'results_table');
       

   catch ME
     warning('Error processing file %d (%s): %s', idx, all_mat_files(idx).name, ME.message);

        % Extract study, participant, and error message
        error_i = results_row.i;  % File index
        error_study = results_row.study;  % Study name
        error_participant = results_row.participant;  % Participant name
        error_message = string(ME.message);  % Error message as string

        % Add the new error row to the table
        new_error = {error_i, error_study, error_participant, error_message};
        error_log = [error_log; new_error];

        % Save the error log after each error, in case the script stops
        save(error_log_file, 'error_log');

        % Continue to the next file in the loop
        continue;
    end
    time = toc;
    fprintf('Duration for this file = %.2f\n', time);

end

% Compute arousals (stage shifts from non-zero to zero divided by the SPT)
% Step 1: Convert SPT to hours
SPT_column = results_table.SPT;  % Assuming 'SPT' is in minutes
SPT_hours = SPT_column / 60;     % Convert SPT to hours
% Step 2: Sum all 'count_X_to_0' columns, excluding 'count_0_to_0'
count_columns = results_table{:, contains(results_table.Properties.VariableNames, 'count_') & ...
                                     contains(results_table.Properties.VariableNames, '_to_0') & ...
                                     ~contains(results_table.Properties.VariableNames, 'count_0_to_0')};
total_count_X_to_0 = sum(count_columns, 2, 'omitnan');  % Sum across columns for each row (participant)
% Step 3: Calculate wakeups per hour
wakeups_per_hour = total_count_X_to_0 ./ SPT_hours;
% Step 4: Append the new column to 'results_table'
results_table.wakeups_per_hour = wakeups_per_hour;
% Save the updated results_table
save(results_table_file, 'results_table');

 % Write to CSV in base_folder
 csv_filename = fullfile(base_folder, 'results_table.csv');
 writetable(results_table, csv_filename);
       