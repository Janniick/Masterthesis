%% add functions
addpath(genpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Feature extraction\Feature_extraction_functions'));
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
% py.importlib.import_module('yasa');
%% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data\THA\';  % Replace with your folder path
% Search for all files ending with '_preprocessed.mat' in all subdirectories
mat_files = dir(fullfile(base_folder, '**', '*_preprocessed_EOG.mat'));

% Initialize the results table with titles
alignment_results = cell(length(mat_files) + 1, 6);  % Adjust to 6 columns

% Set the titles for each column
alignment_results(1, :) = {'Filename', 'EEG Data Duration (s)', 'Scoring Duration (s)', 'Mismatch (s)', 'Mismatch (epochs)', 'Better Alignment'};
%%
for i = 1:length(mat_files)
    % load the file
    load([mat_files(i).folder, '\', mat_files(i).name]);
    % extract the name of the file
    filename = mat_files(i).name;
    fprintf(filename);
    disp(i);
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
    % tolerance = 0.01 * data_duration;  % 1% tolerance
    % if abs(data_duration - scoring_duration) > tolerance
    %     error(['Duration mismatch: The scoring duration is not within 1% of the data duration. Stopping execution.']);
    % end
    % Check for clean epochs
    % Generates EEG_clean.clean_epochs
    % EEG_clean = scoring_clean_epochs(EEG_clean);
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

    % upsample
    % Repeat each entry in hypnogram for 30 * srate times
    epoch_duration_sec = 30;
    repeats_per_epoch = epoch_duration_sec * srate;  % Number of samples per epoch
    hypnogram_upsampled = repelem(hypnogram, repeats_per_epoch);

    % Get the duration of EEG_clean.data and scoring_long_30s in seconds
    data_duration_seconds = size(EEG_clean.data, 2) / EEG_clean.srate;
    scoring_duration_seconds = height(EEG_clean.scoring_long_30s) * 30;  % Assuming scoring_long_30s is in 30s epochs

    % Calculate the mismatch in seconds
    mismatch_seconds = abs(data_duration_seconds - scoring_duration_seconds);
    mismatch_epochs = mismatch_seconds/30;

    % Extract C3 channel data (or use a representative channel for alignment)
    channel_idx_C3 = find(strcmpi({EEG_clean.chanlocs.labels}, 'C3'));
    eeg_data_C3 = EEG_clean.data(channel_idx_C3, :);  % Extract C3 channel data


    % Check if EEG data is shorter than the hypnogram
    if length(eeg_data_C3) < length(hypnogram_upsampled)
        warning('EEG data is shorter than the hypnogram. Skipping alignment.');
        better_alignment = 'N/A';  % You can mark this as 'N/A' or any other value indicating it was skipped
    else
        % Extract the lengths for comparison
        eeg_data_start = eeg_data_C3(1:length(hypnogram_upsampled));  % Align hypnogram at the start
        eeg_data_end = eeg_data_C3(end-length(hypnogram_upsampled)+1:end);  % Align hypnogram at the end

        % Compute correlation for start alignment
        corr_start = corr(eeg_data_start', hypnogram_upsampled');

        % Compute correlation for end alignment
        corr_end = corr(eeg_data_end', hypnogram_upsampled');

        % Determine which alignment is better
        if corr_start >= corr_end
            better_alignment = 'start';
        else
            better_alignment = 'end';
        end
    end


    % Store the results for this file
    alignment_results{i + 1, 1} = mat_files(i).name;  % Filename
    alignment_results{i + 1, 2} = data_duration_seconds;  % Duration of EEG data
    alignment_results{i + 1, 3} = scoring_duration_seconds;  % Duration of scoring data
    alignment_results{i + 1, 4} = mismatch_seconds;  % Mismatch in seconds
    alignment_results{i + 1, 5} = mismatch_epochs;  % Mismatch in epochs
    % Assign the alignment result only if 'better_alignment' is not 'N/A'
    if ~strcmp(better_alignment, 'N/A')
        alignment_results{i + 1, 6} = better_alignment;
    else
        alignment_results{i + 1, 6} = '0';  % Or any placeholder for skipped alignment
    end
    % Better alignment (start/end)
end


% Convert alignment_results to a table before saving
alignment_table = cell2table(alignment_results(2:end, :), 'VariableNames', alignment_results(1, :));

% Create the full file path to save the CSV file in base_folder
csv_file_path = fullfile(base_folder, 'alignment_results.csv');

% Save the table as a CSV file
writetable(alignment_table, csv_file_path);

disp(['Alignment results saved to: ' csv_file_path]);

