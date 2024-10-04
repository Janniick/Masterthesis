%% add functions
addpath(genpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Feature extraction\Feature_extraction_functions'));
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
% py.importlib.import_module('yasa');

%% load python
pyenv('Version', 'C:\Users\janni\AppData\Local\Programs\Python\Python311\python.EXE');
pyrun("import yasa", "yasa");
%% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data';  % Replace with your folder path
% Search for all files ending with '_preprocessed.mat' in all subdirectories
all_mat_files = dir(fullfile(base_folder, '**', '*_preprocessed.mat'));
%% bands
delta_b = [0.5, 4];
theta_b = [4, 8];
alpha_b = [8, 12];
beta_b = [15, 30];
gamma_b = [30, 40];
spindle_b = [12, 16];

for i = 1:1%length(all_mat_files)
    % load the file
    load([all_mat_files(i).folder, '\', all_mat_files(i).name]);
    % extract the name of the file
    filename = all_mat_files(i).name;
    % Split the filename by underscores ('_') and periods ('.')
    split_parts = split(filename, {'_', '.'});
    % Extract the study and participant variables
    study = split_parts{1};
    participant = split_parts{2};
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
    % compare data length and scoring length
    data_duration = length(EEG_clean.data)/EEG_clean.srate/60/60;
    scoring_duration = height(EEG_clean.scoring_long_30s)*30/60/60;
    % Display the duration
    disp(['Duration of EEG_clean.data: ', num2str(data_duration)]);
    disp(['Duration of EEG_clean.scoring_long_30s: ', num2str(scoring_duration)]);
    % Check if scoring_duration is within a 1% margin of data_duration
    tolerance = 0.01 * data_duration;  % 1% tolerance
    if abs(data_duration - scoring_duration) > tolerance
        error(['Duration mismatch: The scoring duration is not within 1% of the data duration. Stopping execution.']);
    end
    % Check for clean epochs
    % Generates EEG_clean.clean_epochs
    EEG_clean = scoring_clean_epochs(EEG_clean);
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



    % Making the hypnogram
    % Extract the sleep stages from the whole night
    sleep_stages_whole = EEG_clean.scoring_long_30s.stage; 
    % Determine the number of epochs (each epoch is 30 seconds)
    num_epochs = length(sleep_stages_whole);
    % Initialize the hypnogram vector to store sleep stages per epoch
    hypnogram = zeros(1, num_epochs);
    % Loop over each epoch and assign the corresponding sleep stage to the hypnogram
    % 0=wake, 1=NREM1, 2=NREM2, 3&4=NREM3, 5=REM, 6<=movement
    for epoch = 1:num_epochs
        stage = sleep_stages_whole(epoch);
        if stage == 3 || stage == 4
            % Combine stage 3 and 4 into NREM3 (using 3 as NREM3 code)
            hypnogram(epoch) = 3;
        elseif stage >= 6
            % Anything above 5 (movement) is set to 6
            hypnogram(epoch) = 6;
        else
            % Use the stage as it is for other sleep stages (0 to 5)
            hypnogram(epoch) = stage;
        end
    end
    % Extraction of whole night features
    % Convert MATLAB hypnogram to Python format using py.numpy.array
    hypnogram_py = py.numpy.array(hypnogram);
    % Call YASA's sleep statistics function
    % TIB (Total Time in Bed)
    % SPT (Sleep Period Time): The time between sleep onset and the final awakening
    % WASO (Wake After Sleep Onset):
    % TST (Total Sleep Time): The sum of time spent in all sleep stages (NREM and REM).
    % N1, N2, N3: The duration (in minutes) spent in different non-REM (NREM) sleep stages:
    % REM: The duration (in minutes) spent in Rapid Eye Movement (REM) sleep.
    % NREM: The total duration of all NREM stages (N1 + N2 + N3), measured in minutes.
    % SOL (Sleep Onset Latency): The amount of time (in minutes) it takes to fall asleep from the moment the person goes to bed.
    % Lat_N1, Lat_N2, Lat_N3, Lat_REM: The latency (time taken) to reach each of the respective stages after sleep onset, measured in minutes:
    % %N1, %N2, %N3, %REM, %NREM: The percentage of total sleep time spent in each sleep stage:
    % SE (Sleep Efficiency): The ratio of Total Sleep Time (TST) to Time in Bed (TIB), expressed as a percentage. It shows how efficiently the person is sleeping.
    % SME (Sleep Maintenance Efficiency): The percentage of the Sleep Period Time (SPT) that the person spends asleep. It's another measure of sleep efficiency but excludes wakefulness before sleep onset.
    sleep_stats = py.yasa.sleep_statistics(hypnogram_py, sf_hyp=1/30);
    sleep_stats = dictionary(sleep_stats);



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
    

    
    % Calculate sleep cycles
    % values = sleep stage
    % start = start epoch
    % length = how many epochs
    sleep_cycles = py.yasa.hypno_find_periods(hypnogram_py, sf_hypno = 1/30, equal_length = false);
    % Extract data from the Python sleep_cycles DataFrame
    % Extract data from the Python sleep_cycles DataFrame
    values = double(sleep_cycles.values);  % Sleep stages
    % Initialize a new array to store the collapsed values
    collapsed_values = [];
    current_stage = values(1, 1);
    current_start = values(1, 2);
    current_length = values(1, 3);
    for i = 2:size(values, 1)
        % Check if the stage is the same as the previous one
        if values(i, 1) == current_stage
            % If the stage is the same, just add to the current length
            current_length = current_length + values(i, 3);
        else
            % If the stage changes, save the current stage, start, and length
            collapsed_values = [collapsed_values; current_stage, current_start, current_length];
            % Update the current stage, start, and length to the new one
            current_stage = values(i, 1);
            current_start = values(i, 2);
            current_length = values(i, 3);
        end
    end
    % Add the last stage to the collapsed array
    collapsed_values = [collapsed_values; current_stage, current_start, current_length];
    % Initialize variables
    sleep_cycles = [];
    start_cycle = 0;
    end_cycle = 0;
    % Loop through collapsed values and find sleep cycles
    for i = 1:size(collapsed_values, 1)
        stage = collapsed_values(i, 1);
        start_epoch = collapsed_values(i, 2);
        length_epoch = collapsed_values(i, 3);
        if stage == 5
            % When REM stage (5) is encountered, mark it as the end of a potential cycle
            if start_cycle > 0
                % Wait until the next NREM stage (2 or 3) to define the end of the cycle
                if i < size(collapsed_values, 1) && (collapsed_values(i + 1, 1) == 2 || collapsed_values(i + 1, 1) == 3)
                    end_cycle = collapsed_values(i + 1, 2)-1;  % Start epoch of next NREM
                    % Calculate the duration of the cycle in epochs
                    duration_in_epochs = end_cycle - start_cycle;
                    % Store the start, end, and duration
                    sleep_cycles = [sleep_cycles; start_cycle, end_cycle, duration_in_epochs];
                    % Reset start_cycle for the next sleep cycle
                    start_cycle = 0;
                end
            end
        elseif stage == 2 || stage == 3
            % If the stage is NREM (2 or 3), mark it as the potential start of a cycle
            if start_cycle == 0
                start_cycle = start_epoch;
            end
        end
    end
    % Convert duration from epochs to minutes (assuming each epoch is 30 seconds)
    duration_in_minutes = (sleep_cycles(:, 3) * 30) / 60;
    % Create a table with appropriate column names
    n_cycles = size(sleep_cycles, 1);  % Number of sleep cycles
    start_epochs = sleep_cycles(:, 1);  % Start epochs of each cycle
    end_epochs = sleep_cycles(:, 2);  % End epochs of each cycle
    duration_minutes = (sleep_cycles(:, 3) * 30) / 60;  % Duration of each cycle in minutes
    % Create the table
    cycle_table = table((1:n_cycles)', start_epochs, end_epochs, duration_minutes, ...
        'VariableNames', {'n_cycle', 'start_epoch', 'end_epoch', 'duration_min'});
    


    % Extract the last sleep stage before wakeup
    last_stage_before_wakeup = extract_last_sleep_stage_before_wakeup(hypnogram);



    % wasserstein interhemispheric spectral measure
    % This function computes the Wasserstein distance between two EEG channels based on their power spectral density (PSD).
    % The PSD is calculated using the Welch method over the specified window and overlap parameters. The function then calculates
    % the cumulative distribution of the PSD in the frequency range of 0.5 to 40 Hz and compares the two channels using the
    % Wasserstein distance metric.
    % done over whole night
    % remove no channels
    remove_channel_which = zeros(1,8);
    duration_window = 4*EEG_clean.srate;
    which_window = hanning(duration_window);
    overlap_window = duration_window*0.5;
    % Saves it automatically with appropriate name, e.g. wasserstein_f3f4
    wasserstein_distance(EEG_clean, 'F3', 'F4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
    wasserstein_distance(EEG_clean, 'C3', 'C4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
    wasserstein_distance(EEG_clean, 'P3', 'P4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
    wasserstein_distance(EEG_clean, 'O1', 'O2', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);



    % Calculate sleep 
    ECG_results = calculate_ECG_by_cycles_and_stages(EEG_clean, cycle_table, collapsed_values);



    %%% Calculations only on clean epochs
    % group the clean epochs by sleep stage
    % Extract necessary data from EEG_clean.clean_epochs
    epochs = EEG_clean.clean_epochs.epoch;
    stages = EEG_clean.clean_epochs.stage;
    start_indices = EEG_clean.clean_epochs.start_index;
    end_indices = EEG_clean.clean_epochs.end_index;
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


   
    % Power spectrum analysis
    % Initialize structs to hold power spectrum results
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



    % FOOOF analysis
    % Channels grouped into FC and PO
    fc_channels = find(contains({EEG_clean.chanlocs.labels}, {'F', 'C'}));  % Find FC channels
    po_channels = find(contains({EEG_clean.chanlocs.labels}, {'P', 'O'}));  % Find PO channels
    % Initialize structs to hold FOOOF results for each stage and overall
    fooof_results = struct('NREM1', [], 'NREM2', [], 'NREM3', [], 'REM', [], 'Overall', []);
    bands = {'delta', 'theta', 'alpha', 'beta', 'spindle', 'gamma'};
    overall_fooof_offset = struct('FC', [], 'PO', []);
    overall_fooof_exponent = struct('FC', [], 'PO', []);
    overall_normalized_absolute_power = struct('FC', [], 'PO', []);
    overall_peak_values = struct('FC', [], 'PO', []);
    stage_fooof_offset = struct('FC', [], 'PO', []);
    stage_fooof_exponent = struct('FC', [], 'PO', []);
    stage_normalized_absolute_power = struct('FC', [], 'PO', []);
    stage_peak_values = struct('FC', [], 'PO', []);
    % Initialize arrays to accumulate FOOOF values across all stages for the "Overall" calculation
    overall_fooof_offset.FC = [];
    overall_fooof_offset.PO = [];
    overall_fooof_exponent.FC = [];
    overall_fooof_exponent.PO = [];
    overall_normalized_absolute_power.FC = [];
    overall_normalized_absolute_power.PO = [];
    overall_peak_values.FC = [];
    overall_peak_values.PO = [];
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
                % Accumulate results for the overall calculation (FC)
                overall_fooof_offset.FC = [overall_fooof_offset.FC; fooof_offset];
                overall_fooof_exponent.FC = [overall_fooof_exponent.FC; fooof_exponent];
                overall_normalized_absolute_power.FC = [overall_normalized_absolute_power.FC; struct2array(normalized_absolute_power)];
                overall_peak_values.FC = [overall_peak_values.FC; struct2array(peak_table)];
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
                % Accumulate results for the overall calculation (PO)
                overall_fooof_offset.PO = [overall_fooof_offset.PO; fooof_offset];
                overall_fooof_exponent.PO = [overall_fooof_exponent.PO; fooof_exponent];
                overall_normalized_absolute_power.PO = [overall_normalized_absolute_power.PO; struct2array(normalized_absolute_power)];
                overall_peak_values.PO = [overall_peak_values.PO; struct2array(peak_table)];
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



    % Extract only clean epoch data for SW and Spindle detection
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
    % Convert the clean EEG data to Python-friendly format (NumPy array)
    data_for_py = cell(1, size(EEG_clean.clean_data, 1));
    for row = 1:size(EEG_clean.clean_data, 1)
        data_for_py{row} = EEG_clean.clean_data(row, :);
    end
    data_for_py = py.numpy.array(data_for_py);
    % Prepare hypnogram data (You need to ensure `hypnodata` corresponds to clean epochs)
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
    % Upsample scoring to match sampling frequency of data
    epoch_duration_sec = 30;  % Duration of each epoch in seconds
    srate = EEG_clean.srate;  % Sampling rate of the EEG data
    % Repeat each entry in hypnogram_clean for 30 * srate times
    repeats_per_epoch = epoch_duration_sec * srate;  % Number of samples per epoch
    hypnogram_clean_upsampled = repelem(hypnogram_clean, repeats_per_epoch);



    % Run the slow wave detection in YASA
    coupling_py = py.dict(pyargs('freq_sp', 12:16, 'p', 0.05, 'time', 1));
    sws_detection_results = py.yasa.sw_detect(data = data_for_py, ...
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
            SW_results.(['Cycle_' num2str(cycle_idx)]).mean = NaN;
            SW_results.(['Cycle_' num2str(cycle_idx)]).std = NaN;
        end
    end
    % Calculate the difference between the first and last cycle
    first_cycle = SW_results.Cycle_1;
    last_cycle = SW_results.(['Cycle_' num2str(height(cycle_table))]);
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








    % Run the spindle detection from YASA
    % Correct the 'thresh' parameter to be a Python dictionary
    thresh_py = py.dict(pyargs('corr', 0.65, 'rel_pow', 0.2, 'rms', 1.5));
    spindles_detection_results = py.yasa.spindles_detect(data = data_for_py, ...
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
        % Convert other numeric columns if necessary
        spindles_detection_info.Start = cellfun(@double, spindles_detection_info.Start);
        spindles_detection_info.Peak = cellfun(@double, spindles_detection_info.Peak);
        spindles_detection_info.End = cellfun(@double, spindles_detection_info.End);
        spindles_detection_info.Duration = cellfun(@double, spindles_detection_info.Duration);
        spindles_detection_info.Amplitude = cellfun(@double, spindles_detection_info.Amplitude);
        spindles_detection_info.RMS = cellfun(@double, spindles_detection_info.RMS);
        spindles_detection_info.AbsPower = cellfun(@double, spindles_detection_info.AbsPower);
        spindles_detection_info.RelPower = cellfun(@double, spindles_detection_info.RelPower);
        spindles_detection_info.Frequency = cellfun(@double, spindles_detection_info.Frequency);
        spindles_detection_info.Oscillations = cellfun(@double, spindles_detection_info.Oscillations);
        spindles_detection_info.Symmetry = cellfun(@double, spindles_detection_info.Symmetry);
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
            Spindle_results.(['Cycle_' num2str(cycle_idx)]).mean = NaN;
            Spindle_results.(['Cycle_' num2str(cycle_idx)]).std = NaN;
        end
    end
    % Calculate the difference between the first and last cycle
    first_cycle = Spindle_results.Cycle_1;
    last_cycle = Spindle_results.(['Cycle_' num2str(height(cycle_table))]);
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






end