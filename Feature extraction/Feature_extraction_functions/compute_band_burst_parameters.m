% Function to compute burst parameters for multiple EEG frequency bands using dual threshold burst detection
% 
% This function computes burst characteristics (number of bursts, mean duration and duration spread) for 
% multiple frequency bands (theta, alpha, beta, spindle) in EEG data. It uses the dual threshold method 
% for burst detection and computes various burst statistics for each frequency band.
% 
% Inputs:
% - EEG_clean: Struct containing EEG data, including the data matrix and sampling rate.
% - chan: The index of the EEG channel to analyze.
% - indexs: The range of data indices to analyze within the EEG data.
% - theta_b: Frequency range for theta bursts (e.g., [4, 8] Hz).
% - alpha_b: Frequency range for alpha bursts (e.g., [8, 13] Hz).
% - beta_b: Frequency range for beta bursts (e.g., [13, 30] Hz).
% - spindle_b: Frequency range for spindle bursts (e.g., [12, 16] Hz).
%
% Outputs:
% - burst_stats: A struct containing the number of bursts, mean burst duration, and duration spread 
%   (standard deviation of burst intervals) for each frequency band. The struct contains the following fields:
%   - burst_band_n: Number of bursts in certain band
%   - burst_band_duration_mean: Mean duration of bursts in certain band
%   - burst_band_duration_ssd: Spread of burst durations in certain band (standard deviation)
%   - burst_band_ssd: Spread of midpoints between burst events in certain band
%    
% Example usage:
% burst_stats = compute_band_burst_parameters(EEG_clean, chan, indexs, theta_b, alpha_b, beta_b, spindle_b);


function burst_stats = compute_band_burst_parameters(EEG_clean, chan, indexs, theta_b, alpha_b, beta_b, spindle_b)

    % Define frequency bands and parameters for burst detection
    freq_bands = {'theta', theta_b; 'alpha', alpha_b; 'beta', beta_b; 'spindle', spindle_b};
    dual_thresh = [1, 1.5];  % Threshold for burst detection across all bands

    % Initialize output structure
    burst_stats = struct();

    % Loop through each frequency band
    for i = 1:size(freq_bands, 1)
        band_name = freq_bands{i, 1};
        band_range = freq_bands{i, 2};

        % Detect bursts using dual threshold method
        burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(EEG_clean.data(chan, indexs)), ...
            fs=EEG_clean.srate, dual_thresh=dual_thresh, f_range=band_range);

        % Compute burst events
        burst_diff = diff([0, double(burst.tolist()), 0]);

        % Compute spread of midpoints of burst events
        if sum(abs(burst_diff)) / 2 > 2
            midpp = (find(burst_diff == 1) + find(burst_diff == -1) - 1) / 2;
            burst_ssd = std(diff(midpp) / EEG_clean.srate);
        else
            burst_ssd = NaN;
        end

        % Compute burst statistics
        burst_stat = py.neurodsp.burst.compute_burst_stats(burst, fs=EEG_clean.srate);
        burst_stat = struct(burst_stat);
        burst_stat.n_bursts = double(burst_stat.n_bursts);
        burst_stat.duration_mean = double(burst_stat.duration_mean);
        burst_stat.duration_std = double(burst_stat.duration_std);

        % Store results into specific variable names for each band
        burst_stats.(sprintf('burst_%s_n', band_name)) = burst_stat.n_bursts;
        burst_stats.(sprintf('burst_%s_duration_mean', band_name)) = burst_stat.duration_mean;
        burst_stats.(sprintf('burst_%s_duration_ssd', band_name)) = burst_stat.duration_std;
        burst_stats.(sprintf('burst_%s_ssd', band_name)) = burst_ssd;  % Spread of midpoints
    end
end
