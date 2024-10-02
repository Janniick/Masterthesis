%% add functions
addpath(genpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Feature extraction\Feature_extraction_functions'));
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
% py.importlib.import_module('yasa');

%% load python
pyenv('Version', 'C:\Users\janni\AppData\Local\Programs\Python\Python311\python.EXE');
%% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data';  % Replace with your folder path
% Search for all files ending with '_preprocessed.mat' in all subdirectories
all_mat_files = dir(fullfile(base_folder, '**', '*_preprocessed.mat'));
%% load the file
i=38;
load([all_mat_files(i).folder, '\', all_mat_files(i).name]);
%% extract the name of the file
% Example filename
filename = all_mat_files(i).name;
% Split the filename by underscores ('_') and periods ('.')
split_parts = split(filename, {'_', '.'});
% Extract the study and participant variables
study = split_parts{1};      
participant = split_parts{2}; 
%% Expand the scoring if it's condensed
% Generates EEG_clean.scoring_long
EEG_clean = scoring_expansion(EEG_clean);
%% Standardize the scoring
% Only keep epoch, stage and duration
% Generates EEG_clean.scoring_long_standardized
EEG_clean = scoring_standardization(EEG_clean);
%% convert 20s to 30s duration if needed
% Generates EEG_clean.scoring_long_30s
EEG_clean = scoring_convert20sTo30s(EEG_clean);
%% add indexes of every epoch to the scoring
% Updates EEG_clean.scoring_long_30s
EEG_clean = scoring_add_start_end_indices(EEG_clean);
%% compare data length and scoring length
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
%% Check for clean epochs
% Generates EEG_clean.clean_epochs
EEG_clean = scoring_clean_epochs(EEG_clean);
%% initialize things for spectrum and so on
% stfft sepctrogram
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

%% bands
delta_b = [0.5, 4];
theta_b = [4, 8];
alpha_b = [8, 12];
beta_b = [15, 30];
gamma_b = [30, 40];
spindle_b = [12, 16];

%% loop over channels
nr_chan = length(EEG_clean.chanlocs);
% tic;  % Start timing
nr_epochs = height(EEG_clean.clean_epochs);

% Preallocate a cell array for each worker's results
parfor_results = cell(nr_chan, 1);  % One cell per channel
% Open a parallel pool if it isn't already open
if isempty(gcp('nocreate'))
    parpool(max(1, floor(0.8 * feature('numcores'))));  % This opens a parallel pool with 0.8*available workers
end
parfor chan = 1:1%nr_chan %loops over all channels
    current_channel = EEG_clean.chanlocs(chan).labels;
    temp_table = cell(nr_epochs, 1); % Temporary storage for each channel within this parfor loop
    %loop over all clean epochs
    for epochs = 1:10%nr_epochs
        current_epoch = EEG_clean.clean_epochs.epoch(epochs);
        start_index = EEG_clean.clean_epochs{epochs, 'start_index'}; %(epochs-1)*epoch_windowlen*(1-epoch_windowover)*EEG_clean.srate + start_index_o; % (epochs-1)*epoch_windowlen*EEG_clean.srate + start_index_o;
        end_index = EEG_clean.clean_epochs{epochs, 'end_index'}; % min(start_index + epoch_windowlen*EEG_clean.srate, end_index_o) % min(epochs*epoch_windowlen*EEG_clean.srate + start_index_o, end_index_o);
        indexs = start_index:end_index;
        % Check if the 'stage' is a string, convert to a number if true, otherwise keep it as is
        if ischar(EEG_clean.clean_epochs.stage(epochs)) || isstring(EEG_clean.clean_epochs.stage(epochs))
            current_stage = str2num(EEG_clean.clean_epochs.stage(epochs));
        else
            current_stage = EEG_clean.clean_epochs.stage(epochs);  % Leave it as is if it's already numeric
        end

        % Power Spectrum: pwelch
        % Calculate the absolute and relative power of certain frequency bands using the pwelch method
        tic;
        [absolute_power, relative_power, freq, psd] = calculate_power_spectrum(EEG_clean, chan, indexs, srate, freq_res, delta_b, theta_b, alpha_b, beta_b, spindle_b);
        duration_spectrum = toc;
        disp(duration_spectrum);
        % fooof algorithm (Fitting Oscillations and One-Over-F)
        % analyze neural power spectra by separating oscillatory components from the aperiodic background activity
        tic;
        [fooof_offset, fooof_exponent, peak_table, normalized_absolute_power] = run_fooof_and_detect_strongest_peaks(freq, psd, delta_b, theta_b, alpha_b, beta_b, spindle_b, gamma_b);
        duration_fooof = toc;
        disp(duration_fooof);

        % phase amp coupling
        % a measure of the interaction between the phase of low-frequency brain oscillations and the amplitude of higher-frequency oscillations
        % PACTool Matlab Add on
        % https://github.com/sccn/PACTools/tree/develop
        tic;
        % Create a new EEG structure with only the necessary fields
        EEG_cleant = struct();  % Create a new empty structure for EEG_cleant
        EEG_cleant.srate = EEG_clean.srate;  % Copy over only the required fields
        EEG_cleant.chanlocs = EEG_clean.chanlocs;  % Copy channel locations
        EEG_cleant.data = EEG_clean.data(:, indexs);  % Copy the sliced data
        EEG_cleant.times = EEG_clean.times(:, indexs);  % Copy the sliced times
        EEG_cleant.etc = EEG_clean.etc;
        EEG_cleant.pnts = length(indexs);  % Number of points in the current epoch
        EEG_cleant.trials = 1;  % Single epoch
        EEG_cleant = pop_pac_modified(EEG_cleant, 'channels', delta_b, spindle_b, chan, chan, 'forcecomp', 1);
        phase_amplitude_coupling_delta_spindle = mean(EEG_cleant.etc.eegpac(1).glm.pacval, 'all');
        EEG_cleant.etc.eegpac = [];
        EEG_cleant.etc.pacplotopt = [];

        duration_coupl = toc;
        disp(duration_coupl);

        % Complexity measures
        % Entropy calculations
        % https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20-%20MatLab/README.md
        tic;
        perm_en = PermEn(symb_king(EEG_clean.data(chan, indexs)), 7);
        perm_en = perm_en(7);
        perm_en_t = PermEn(symb_king(theta_ts), 7);
        perm_en_t = perm_en_t(7);
        perm_en_a = PermEn(symb_king(alpha_ts), 7);
        perm_en_a = perm_en_a(7);
        perm_en_b = PermEn(symb_king(beta_ts), 7);
        perm_en_b = perm_en_b(7);
        perm_en_s = PermEn(symb_king(spindle_ts), 7);
        perm_en_s = perm_en_s(7);

        % hilbert envelope hurst etc
        hil = hilbert(EEG_clean.data(chan, indexs));
        hil_env = abs(hil);
        hurst_hilbert = hurst_estimate(hil_env, 'absval', 0, 1);
        slope_en_hilbert =  SlopEn(zscore(hil_env));
        std_hilbert =  std(hil_env);
        mean_hilbert =  mean(hil_env);
        rmssd_hilbert = sqrt(mean(diff(hil_env).^2));
        lzw_hilbert = lzw(hil_env>mean_hilbert);

        % Lagged coherence is a measure to quantify the rhythmicity of neural signals.
        % Lagged coherence works by quantifying phase consistency between non-overlapping data fragments, calculated with Fourier coefficients.
        % The consistency of the phase differences across epochs indexes the rhythmicity of the signal.
        % https://neurodsp-tools.github.io/neurodsp/auto_tutorials/rhythm/plot_LaggedCoherence.html

        t_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=theta_b);
        a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=alpha_b);
        t_a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=[theta_b(1), alpha_b(2)]);
        b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=beta_b);
        a_b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=[alpha_b(1), beta_b(2)]);
        s_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=spindle_b);
        d_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=delta_b);
        d_s_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=[delta_b(1), spindle_b(2)]);

        % The Hurst exponent (or Hurst estimate) is a measure used to evaluate the long-term memory or self-similarity of a time series.
        % It helps you understand the persistence, randomness, or mean-reverting behavior of a dataset.
        % H < 0.5 (Anti-persistent behavior): The time series has anti-persistent or mean-reverting behavior. If the time series increases in the past, it is likely to decrease in the future and vice versa.
        % H = 0.5 (Random walk): The time series is a random walk (Brownian motion), meaning future values are independent of past values. There is no long-term memory in the process.
        % H > 0.5 (Persistent behavior): The time series shows persistent or trending behavior. If it increases in the past, it is likely to keep increasing in the future (and similarly for decreases).
        hurs_exp_long = hurst_estimate(zscore(EEG_clean.data(chan, indexs)), 'absval', 0, 1);
        hurs_alpha = hurst_estimate(zscore(alpha_ts), 'absval', 0, 1);
        hurs_beta = hurst_estimate(zscore(beta_ts), 'absval', 0, 1);
        hurs_theta = hurst_estimate(zscore(theta_ts), 'absval', 0, 1);
        hurs_delta = hurst_estimate(zscore(delta_ts), 'absval', 0, 1);
        hurs_gamma = hurst_estimate(zscore(gamma_ts), 'absval', 0, 1);
        hurs_spindle = hurst_estimate(zscore(spindle_ts), 'absval', 0, 1);

        % Complexity measures
        higuchi_fd = Higuchi_FD(zscore(EEG_clean.data(chan, indexs)), min(50, max(1,length(indexs)-2)));
        katz_fd = Katz_FD(zscore(EEG_clean.data(chan, indexs)));
        bubble_en = BubbEn(zscore(EEG_clean.data(chan, indexs)));
        spect_en = SpecEn(zscore(EEG_clean.data(chan, indexs)));
        sydy_en = SyDyEn(zscore(EEG_clean.data(chan, indexs)));
        phase_en = PhasEn(zscore(EEG_clean.data(chan, indexs)));
        slope_en =  SlopEn(EEG_clean.data(chan, indexs)); % original data
        slope_en_z =  SlopEn(zscore(EEG_clean.data(chan, indexs))); %zscore
        slope_en_z_05_20 =  SlopEn(zscore(EEG_clean.data(chan, indexs)), 'Lvls',[0.5,20]); %zscore
        eoe_en = EnofEn(zscore(EEG_clean.data(chan, indexs)));
        att_en = AttnEn(zscore(EEG_clean.data(chan, indexs)));

        [MSx,smse] = MSEn(EEG_clean.data(chan, indexs), MSobject('SampEn')); % multiscale sample enntropy

        scalz = unique(2*floor(linspace(EEG_clean.srate/30, EEG_clean.srate/0.7, 5)/2)+ 1);
        mlz = getMultiscaleLZ(EEG_clean.data(chan, indexs), scalz);
        smlz = sum(mlz);

        % other time vars
        kurt = kurtosis(EEG_clean.data(chan, indexs));
        skew = skewness(EEG_clean.data(chan, indexs));
        [hjorth_activity, hjorth_mobility, hjorth_complexity] = hjorth(EEG_clean.data(chan, indexs)',0);

        duration_comp = toc;
        disp(duration_comp);

        %  save the data jannick
        % Add Burst Features
        burst_theta_n = burst_parameters.burst_theta_n;
        burst_theta_duration_mean = burst_parameters.burst_theta_duration_mean;
        burst_theta_duration_ssd = burst_parameters.burst_theta_duration_ssd;
        burst_alpha_n = burst_parameters.burst_alpha_n;
        burst_alpha_duration_mean = burst_parameters.burst_alpha_duration_mean;
        burst_alpha_duration_ssd = burst_parameters.burst_alpha_duration_ssd;
        burst_beta_n = burst_parameters.burst_beta_n;
        burst_beta_duration_mean = burst_parameters.burst_beta_duration_mean;
        burst_beta_duration_ssd = burst_parameters.burst_beta_duration_ssd;
        burst_spindle_n = burst_parameters.burst_spindle_n;
        burst_spindle_duration_mean = burst_parameters.burst_spindle_duration_mean;
        burst_spindle_duration_ssd = burst_parameters.burst_spindle_duration_ssd;
        % Add Power Features
        absolute_delta = absolute_power.delta;
        absolute_theta = absolute_power.theta;
        absolute_alpha = absolute_power.alpha;
        absolute_beta = absolute_power.beta;
        absolute_spindle = absolute_power.spindle;
        relative_delta = relative_power.delta;
        relative_theta = relative_power.theta;
        relative_alpha = relative_power.alpha;
        relative_beta = relative_power.beta;
        relative_spindle = relative_power.spindle;
        % Add Fooof Features
        fooof_offset_val = fooof_offset;
        fooof_exponent_val = fooof_exponent;
        % Add Peak Table Features
        % Extract values from the structured peak_table
        peak_delta_frequency = peak_table.delta_frequency;
        peak_delta_amplitude = peak_table.delta_amplitude;
        peak_delta_bandwidth = peak_table.delta_bandwidth;
        peak_theta_frequency = peak_table.theta_frequency;
        peak_theta_amplitude = peak_table.theta_amplitude;
        peak_theta_bandwidth = peak_table.theta_bandwidth;
        peak_alpha_frequency = peak_table.alpha_frequency;
        peak_alpha_amplitude = peak_table.alpha_amplitude;
        peak_alpha_bandwidth = peak_table.alpha_bandwidth;
        peak_beta_frequency = peak_table.beta_frequency;
        peak_beta_amplitude = peak_table.beta_amplitude;
        peak_beta_bandwidth = peak_table.beta_bandwidth;
        peak_spindle_frequency = peak_table.spindle_frequency;
        peak_spindle_amplitude = peak_table.spindle_amplitude;
        peak_spindle_bandwidth = peak_table.spindle_bandwidth;
        peak_gamma_frequency = peak_table.gamma_frequency;
        peak_gamma_amplitude = peak_table.gamma_amplitude;
        peak_gamma_bandwidth = peak_table.gamma_bandwidth;
        % Add Normalized Absolute Power Features
        norm_absolute_delta = normalized_absolute_power.delta;
        norm_absolute_theta = normalized_absolute_power.theta;
        norm_absolute_alpha = normalized_absolute_power.alpha;
        norm_absolute_beta = normalized_absolute_power.beta;
        norm_absolute_spindle = normalized_absolute_power.spindle;
        norm_absolute_gamma = normalized_absolute_power.gamma;
        norm_absolute_total = normalized_absolute_power.total;
        % Add Mutual Information Band Features
        mi_delta_theta = mutual_information_band.delta_theta;
        mi_delta_alpha = mutual_information_band.delta_alpha;
        mi_delta_beta = mutual_information_band.delta_beta;
        mi_delta_spindle = mutual_information_band.delta_spindles;
        mi_theta_alpha = mutual_information_band.theta_alpha;
        mi_theta_beta = mutual_information_band.theta_beta;
        mi_theta_spindle = mutual_information_band.theta_spindles;
        mi_alpha_beta = mutual_information_band.alpha_beta;
        mi_alpha_spindle = mutual_information_band.alpha_spindles;
        mi_beta_spindle = mutual_information_band.beta_spindles;
        % Add Phase Amplitude Coupling Features
        pac_theta_alpha = phase_amplitude_coupling_theta_alpha;
        pac_alpha_beta = phase_amplitude_coupling_alpha_beta;
        pac_theta_beta = phase_amplitude_coupling_theta_beta;
        pac_delta_spindle = phase_amplitude_coupling_delta_spindle;
        % Add Complexity Measures
        complexity_perm_en = perm_en;
        complexity_hurst = hurs_exp_long;
        complexity_hurst_alpha = hurs_alpha;
        complexity_hurst_beta = hurs_beta;
        complexity_hurst_theta = hurs_theta;
        complexity_hurst_delta = hurs_delta;
        complexity_hurst_gamma = hurs_gamma;
        complexity_hurst_spindle = hurs_spindle;
        complexity_hilbert_rmssd = rmssd_hilbert;
        % Add Lagged Coherence Features
        lag_coh_theta = t_lagcoh;
        lag_coh_alpha = a_lagcoh;
        lag_coh_theta_alpha = t_a_lagcoh;
        lag_coh_beta = b_lagcoh;
        lag_coh_alpha_beta = a_b_lagcoh;
        lag_coh_spindle = s_lagcoh;
        lag_coh_delta = d_lagcoh;
        lag_coh_delta_spindle = d_s_lagcoh;
        % Add Other Features
        kurt_val = kurt;
        skew_val = skew;
        hjorth_activity_val = hjorth_activity;
        hjorth_mobility_val = hjorth_mobility;
        hjorth_complexity_val = hjorth_complexity;
        % Ensure text fields are wrapped in cell arrays or use string arrays for text data
        % Create a new row with all the extracted features
        new_row = table({study}, {participant}, {current_channel}, current_epoch, start_index, end_index, current_stage, ...
            burst_theta_n, burst_theta_duration_mean, burst_theta_duration_ssd, ...
            burst_alpha_n, burst_alpha_duration_mean, burst_alpha_duration_ssd, ...
            burst_beta_n, burst_beta_duration_mean, burst_beta_duration_ssd, ...
            burst_spindle_n, burst_spindle_duration_mean, burst_spindle_duration_ssd, ...
            absolute_delta, absolute_theta, absolute_alpha, absolute_beta, absolute_spindle, ...
            relative_delta, relative_theta, relative_alpha, relative_beta, relative_spindle, ...
            fooof_offset_val, fooof_exponent_val, ...
            peak_delta_frequency, peak_delta_amplitude, peak_delta_bandwidth, ...
            peak_theta_frequency, peak_theta_amplitude, peak_theta_bandwidth, ...
            peak_alpha_frequency, peak_alpha_amplitude, peak_alpha_bandwidth, ...
            peak_beta_frequency, peak_beta_amplitude, peak_beta_bandwidth, ...
            peak_spindle_frequency, peak_spindle_amplitude, peak_spindle_bandwidth, ...
            peak_gamma_frequency, peak_gamma_amplitude, peak_gamma_bandwidth, ...
            norm_absolute_delta, norm_absolute_theta, norm_absolute_alpha, norm_absolute_beta, ...
            norm_absolute_spindle, norm_absolute_gamma, norm_absolute_total, ...
            mi_delta_theta, mi_delta_alpha, mi_delta_beta, mi_delta_spindle, ...
            mi_theta_alpha, mi_theta_beta, mi_theta_spindle, ...
            mi_alpha_beta, mi_alpha_spindle, mi_beta_spindle, ...
            pac_theta_alpha, pac_alpha_beta, pac_theta_beta, pac_delta_spindle, ...
            complexity_perm_en, complexity_hurst, complexity_hurst_alpha, complexity_hurst_beta, ...
            complexity_hurst_theta, complexity_hurst_delta, complexity_hurst_gamma, complexity_hurst_spindle, ...
            complexity_hilbert_rmssd, ...
            lag_coh_theta, lag_coh_alpha, lag_coh_theta_alpha, lag_coh_beta, lag_coh_alpha_beta, ...
            lag_coh_spindle, lag_coh_delta, lag_coh_delta_spindle, ...
            kurt_val, skew_val, hjorth_activity_val, hjorth_mobility_val, hjorth_complexity_val, ...
            'VariableNames', {'study', 'participant', 'channel', 'epoch', 'start_index', 'end_index', 'stage', ...
            'burst_theta_n', 'burst_theta_duration_mean', 'burst_theta_duration_ssd', ...
            'burst_alpha_n', 'burst_alpha_duration_mean', 'burst_alpha_duration_ssd', ...
            'burst_beta_n', 'burst_beta_duration_mean', 'burst_beta_duration_ssd', ...
            'burst_spindle_n', 'burst_spindle_duration_mean', 'burst_spindle_duration_ssd', ...
            'absolute_delta_power', 'absolute_theta_power', 'absolute_alpha_power', 'absolute_beta_power', 'absolute_spindle_power', ...
            'relative_delta_power', 'relative_theta_power', 'relative_alpha_power', 'relative_beta_power', 'relative_spindle_power', ...
            'fooof_offset', 'fooof_exponent', ...
            'peak_delta_frequency', 'peak_delta_amplitude', 'peak_delta_bandwidth', ...
            'peak_theta_frequency', 'peak_theta_amplitude', 'peak_theta_bandwidth', ...
            'peak_alpha_frequency', 'peak_alpha_amplitude', 'peak_alpha_bandwidth', ...
            'peak_beta_frequency', 'peak_beta_amplitude', 'peak_beta_bandwidth', ...
            'peak_spindle_frequency', 'peak_spindle_amplitude', 'peak_spindle_bandwidth', ...
            'peak_gamma_frequency', 'peak_gamma_amplitude', 'peak_gamma_bandwidth', ...
            'normalized_delta_power', 'normalized_theta_power', 'normalized_alpha_power', 'normalized_beta_power', ...
            'normalized_spindle_power', 'normalized_gamma_power', 'normalized_total_power', ...
            'mi_delta_theta', 'mi_delta_alpha', 'mi_delta_beta', 'mi_delta_spindles', 'mi_theta_alpha', 'mi_theta_beta', ...
            'mi_theta_spindles', 'mi_alpha_beta', 'mi_alpha_spindles', 'mi_beta_spindles', ...
            'phase_amplitude_coupling_theta_alpha', 'phase_amplitude_coupling_alpha_beta', ...
            'phase_amplitude_coupling_theta_beta', 'phase_amplitude_coupling_delta_spindle', ...
            'complexity_perm_en', 'complexity_hurst', 'complexity_hurst_alpha', 'complexity_hurst_beta', ...
            'complexity_hurst_theta', 'complexity_hurst_delta', 'complexity_hurst_gamma', 'complexity_hurst_spindle', ...
            'complexity_hilbert_rmssd', ...
            'lagged_coherence_theta', 'lagged_coherence_alpha', 'lagged_coherence_theta_alpha', ...
            'lagged_coherence_beta', 'lagged_coherence_alpha_beta', 'lagged_coherence_spindle', ...
            'lagged_coherence_delta', 'lagged_coherence_delta_spindle', ...
            'kurtosis', 'skewness', 'hjorth_activity', 'hjorth_mobility', 'hjorth_complexity'});

        % Store the row result in the temporary array
        temp_table{epochs} = new_row;

        duration_save = toc;
        disp(duration_save);
    end % epoch loop
    % After processing all epochs, store this channel's table in the worker's result cell
    parfor_results{chan} = vertcat(temp_table{:});
end % channel loop
% After the parfor loop, concatenate the results from all workers
combined_table = vertcat(parfor_results{:});

% Generate the file name dynamically using study and participant variables
file_name = sprintf('%s_%s_features_chan_epoch.csv', study, participant);
% Create the full file path
full_file_path = fullfile(pathtosave, file_name);
% Save the combined table as a CSV file with the dynamic name in the specified directory
writetable(combined_table, full_file_path);
% Display a success message
disp(['File saved successfully: ', full_file_path]);

% elapsedTimeInSeconds = toc;  % End timing and get elapsed time in seconds
% Convert to minutes
elapsedTimeInMinutes = elapsedTimeInSeconds / 60;
% Display the elapsed time in both seconds and minutes
fprintf('The for loop took %.4f seconds (%.4f minutes) to execute.\n', elapsedTimeInSeconds, elapsedTimeInMinutes);

delete(gcp('nocreate'));% This will shut down the parallel pool if it's running


        %% wasserstein interhemispheric spectral measure
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

        wasserstein_distance(EEG_clean, 'F3', 'F4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'C3', 'C4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'P3', 'P4', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);
        wasserstein_distance(EEG_clean, 'O1', 'O2', remove_channel_which, which_window, overlap_window, duration_window, EEG_clean.srate);


        %% ECG measures
        ECG_values = calculate_ECG_values(EEG_clean, indexs);

        %% median mutual information between F, C, P & O channels (was bedeuten die Werte?)
        % Define regions and channels
        regions = {'F', 'C', 'O', 'P'};
        channels = {{'F3', 'F4'}, {'C3', 'C4'}, {'O1', 'O2'}, {'P3', 'P4'}};
        % remove no channels
        remove_channel_which = zeros(1,8);

        % Call the function to calculate mutual information topographically
        mutual_information_topo = mutual_information_topographically(EEG_clean, remove_channel_which, indexs, regions, channels);
    end
    %% Macrofeatures
    epoch_windowlen = 30; % seconds
    epoch_windowover = 0; % no overlap
    if length(EEG_clean.data(1, :))/EEG_clean.srate >= 2
        % add stage to each data sample
        current_stages = str2double(string(EEG_clean.scoring_long_30s.stage)');
        lastvalue = find(~isnan(table2array(EEG_clean.scoring_long_30s(:,'duration_s')))');
        lastvalue = lastvalue(end);
        % hypnodata = [repelem(current_stages(1),table2array(current_scoring(1,'start_s'))*EEG_clean.srate), ...
        %             repelem(current_stages(1:lastvalue), 30*EEG_clean.srate), ...
        %             repelem(current_stages(lastvalue), max(0,(size(EEG_clean.data,2)/EEG_clean.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG_clean.srate))
        %             ];
        % current_epoch = str2double(string(current_scoring.epoch)');
        % epochdata = [repelem(current_epoch(1),table2array(current_scoring(1,'start_s'))*EEG_clean.srate), ...
        %             repelem(current_epoch(1:lastvalue), 30*EEG_clean.srate), ...
        %             repelem(current_epoch(lastvalue), max(0,(size(EEG_clean.data,2)/EEG_clean.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG_clean.srate))
        %             ];
        % Initialize an empty array for hypnodata
        hypnodata = [];

        % Loop through each row of the EEG_clean.scoring_long_30s table
        for i = 1:height(EEG_clean.scoring_long)
            % Get the number of samples for this epoch
            num_samples_in_epoch = EEG_clean.scoring_long.duration_s(1) *srate;

            % Get the stage for this epoch (converted to a number if necessary)
            stage = str2double(EEG_clean.scoring_long.stage{i});

            % Repeat the stage value for the number of samples in the epoch
            hypnodata = [hypnodata, repelem(stage, num_samples_in_epoch)];
        end

        % Check if the lengths of hypnodata and EEG_clean.data match
        data_length = size(EEG_clean.data, 2);  % Number of samples in EEG_clean.data
        hypnodata_length = length(hypnodata);   % Length of hypnodata

        % If hypnodata is shorter than EEG_clean.data
        if hypnodata_length < data_length
            % Calculate the difference in length
            length_diff = data_length - hypnodata_length;

            % Shorten EEG_clean.data by removing the extra samples from the beginning
            EEG_clean.data_shortened = EEG_clean.data(:, length_diff+1:end);

            % Display the result
            disp(['Data was shortened by ', num2str(length_diff), ' samples from the beginning.']);
        else
            EEG_clean.data_shortened = EEG_clean.data;
            disp('hypnodata and EEG_clean.data have the same length. No changes made.');
        end







        data_for_py = cell(1, size(EEG_clean.data_shortened, 1));
        for row = 1:size(EEG_clean.data, 1)
            data_for_py(row) = {EEG_clean.data_shortened(row, :)};
        end
        data_for_py = py.numpy.array(data_for_py);


        sws_detection_result = py.yasa.sw_detect(data = data_for_py, ...
            sf = py.float(EEG_clean.srate), ...
            ch_names=py.list({EEG_clean.chanlocs.labels}), ...
            hypno = py.numpy.array(hypnodata, 'int'), ...
            include = py.numpy.array([2, 3], 'int') , ... % stages to include in SWS detection
            freq_sw = py.tuple([0.5, 4.5]), ...
            dur_neg = py.tuple([0.1, 1.5]), ...
            dur_pos = py.tuple([0.1, 1.5]), ... % duration in seconds
            amp_neg = py.tuple([5, 1000]), ...
            amp_pos = py.tuple([5, 1000]), ...
            amp_ptp = py.tuple([20, 1000]), ...
            coupling = false, ... % phase coupling calculation
            remove_outliers = false, ... % isolation forest to remove wierd detected peaks
            verbose=false);
        py.yasa.plot_hypnogram(hypnodata)

        if ~isempty(sws_detection_info)
            sws_detection_info = sws_detection_result.summary();
            sws_colnames = cellfun(@string,cell(sws_detection_info.columns.values.tolist()));
            sws_detection_info = cell(sws_detection_info.values.tolist());
            sws_detection_info = cellfun(@(xs)cell2table(cell(xs), 'VariableNames', sws_colnames),sws_detection_info, 'UniformOutput', false);
            sws_detection_info = vertcat(sws_detection_info{:});
            sws_detection_info.Stage = cellfun(@double, sws_detection_info.Stage);
            sws_detection_info.Channel = string(sws_detection_info.Channel);
            sws_detection_info.IdxChannel = cellfun(@double, sws_detection_info.IdxChannel);

            sws_detection_info.neg_sw_peaks_ind = sws_detection_info.NegPeak*EEG_clean.srate;
            sws_detection_info.duration_negative_phase = (sws_detection_info.MidCrossing - sws_detection_info.Start)*EEG_clean.srate;
            sws_detection_info.duration_positive_phase = (sws_detection_info.End - sws_detection_info.MidCrossing)*EEG_clean.srate;
            sws_detection_info.duration_start_to_neg = (sws_detection_info.NegPeak - sws_detection_info.Start)*EEG_clean.srate;
            sws_detection_info.duration_neg_to_mid = (sws_detection_info.MidCrossing - sws_detection_info.NegPeak)*EEG_clean.srate;
            sws_detection_info.duration_mid_to_pos = (sws_detection_info.PosPeak - sws_detection_info.MidCrossing)*EEG_clean.srate;
            sws_detection_info.duration_pos_to_end = (sws_detection_info.End - sws_detection_info.PosPeak)*EEG_clean.srate;


            % dbscan
            [npeaks_s, idx_npeaks] = sort(sws_detection_info.NegPeak); % event times as column and in ms (or change), maybe look at ordering as well
            eps = round(1000*(1/4.5+0.15/(0.8))/2)/1000; % max distance in s between two events in the same cluster
            % 1000/4.5 max possible if up to 4.5 Hz,
            % 1000*0.15/(0.8) minimum possible if lowest travelling speed is 0.8 m/s (previously 1.2-7 m/s found) and max 15 cm electrode distance (potentially one or two not capturing the wave)
            minPts = 2; % Minimum number of points to form a dense region
            sw_clusters = dbscan(npeaks_s, eps, minPts); % 2*eps?

            sws_detection_info = sws_detection_info(idx_npeaks,:); % rearrange, so that it corresponds to the same
            sws_detection_info.sw_cluster = sw_clusters;

            % save table per file
            sws_detection_info.id = repelem(string(thenamescore{2}),size(sws_detection_info,1))';
            sws_detection_info.night = repelem(string(thenamescore{3}),size(sws_detection_info,1))';


            % reference peak
            swgr_idx = 1:numel(sw_clusters);

            % representative member
            swgr_i = zeros(1,numel(sw_clusters));
            idx_refwave = floor(unique(accumarray(findgroups(sw_clusters), swgr_idx, [], @mean)));

            if ~isempty(idx_refwave)
                swgr_i(idx_refwave) = 1:length(idx_refwave);
            end
            sws_detection_info.sw_ref_wave = swgr_i';


            sws_detection_info.percent_involved = zeros(size(sws_detection_info,1),1);
            sws_detection_info.mean_slope = zeros(size(sws_detection_info,1),1);
            sws_detection_info.sync_score_bernardi = zeros(size(sws_detection_info,1),1);
            if ~isempty(idx_refwave)
                for refw = 1:length(idx_refwave)
                    ref_neg_idx = sws_detection_info.neg_sw_peaks_ind(idx_refwave(refw));
                    min_20ms = ref_neg_idx - floor(0.02*EEG_clean.srate);
                    plus_20ms = ref_neg_idx + floor(0.02*EEG_clean.srate);

                    sws_detection_info.percent_involved(idx_refwave(refw)) = sum(mean(EEG_clean.data(:, min_20ms:plus_20ms),2)< -5)/size(EEG_clean.data,1);

                    slope1 = abs(sws_detection_info.ValNegPeak(idx_refwave(refw)))/(sws_detection_info.NegPeak(idx_refwave(refw))-sws_detection_info.Start(idx_refwave(refw)));
                    slope2 = abs(sws_detection_info.ValNegPeak(idx_refwave(refw)))/(sws_detection_info.MidCrossing(idx_refwave(refw))-sws_detection_info.NegPeak(idx_refwave(refw)));

                    sws_detection_info.mean_slope(idx_refwave(refw)) = mean([slope1,slope2]);
                    sws_detection_info.sync_score_bernardi(idx_refwave(refw)) = sws_detection_info.percent_involved(idx_refwave(refw))*sws_detection_info.mean_slope(idx_refwave(refw));
                end
            end


            sws_detection_info.name = repelem({thename},size(sws_detection_info,1))';
            sws_detection_info.epoch = epochdata(sws_detection_info.neg_sw_peaks_ind)';
        else
            sws_detection_info = [];
        end



        t_tables{i} =  sws_detection_info;
    end
else
    % save table per file
    t_tables{i} = [];
end





%% Output Variabeln
% power per band
delta_absolute_power
delta_relative_power
theta_absolute_power
theta_relative_power
alpha_absolute_power
alpha_relative_power
beta_absolute_power
beta_relative_power
spindle_absolute_power
spindle_relative_power,
normalized_absolute_delta_power,
normalized_absolute_theta_power,
normalized_absolute_alpha_power,
normalized_absolute_beta_power,
normalized_absolute_gamma_power,
normalized_absolute_spindle_power,
normalized_absolute_total_power, ...
    % band bursts
theta_burst_stat.n_bursts, ...
    theta_burst_stat.duration_mean, ...
    theta_burst_stat.duration_std, ...
    alpha_burst_stat.n_bursts, ...
    alpha_burst_stat.duration_mean, ...
    alpha_burst_stat.duration_std, ...
    beta_burst_stat.n_bursts, ...
    beta_burst_stat.duration_mean, ...
    beta_burst_stat.duration_std, ...
    spindle_burst_stat.n_bursts, ...
    spindle_burst_stat.duration_mean, ...
    spindle_burst_stat.duration_std, ...
    % fooof
fooof_exponent
fooof_offset
strongest_frequency_peak_overall, strongest_frequency_peak_overall_area, ...
    strongest_frequency_peak_delta, strongest_frequency_peak_delta_area, ...
    strongest_frequency_peak_theta, strongest_frequency_peak_theta_area, ...
    strongest_frequency_peak_alpha, strongest_frequency_peak_alpha_area, ...
    strongest_frequency_peak_beta, strongest_frequency_peak_beta_area, ...
    strongest_frequency_peak_spindle, strongest_frequency_peak_spindle_area, ...
    % phase amplitude coupling
phase_amplitude_coupling_theta_alpha, ...
    phase_amplitude_coupling_alpha_beta, ...
    phase_amplitude_coupling_theta_beta, ...
    phase_amplitude_coupling_delta_spindle, ...
    % mutual information
mutual_information_band
mutual_information_topo







%% average pro sleep stage?
%erster fooof löschen weil zweiter detaillierter?
% ideen für sleep EEG?
%yasa spindles?
% cap pattern (cyclic, acyclic pattern)
% arousals über ganze nacht, nicht nur clean epochs?

