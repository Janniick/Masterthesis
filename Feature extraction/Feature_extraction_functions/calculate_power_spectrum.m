% CALCULATE_POWER_SPECTRUM - Computes the absolute and relative power of EEG signals within specified frequency bands.
%
% This function uses the pwelch method to calculate the power spectral density (PSD) of an EEG signal.
% It then calculates the absolute and relative power for the specified frequency bands: delta, theta, alpha, beta, and spindle.
% The frequency and PSD values for the signal are also returned for further analysis.
%
% Inputs:
%   EEG_clean   - EEG dataset structure containing the EEG data
%   chan        - Index of the channel to analyze
%   indexs      - Indices of the EEG data segment to analyze
%   srate       - Sampling rate of the EEG data
%   freq_res    - Frequency resolution for the power spectrum analysis
%   delta_b     - Frequency range for the delta band [start_freq, end_freq]
%   theta_b     - Frequency range for the theta band [start_freq, end_freq]
%   alpha_b     - Frequency range for the alpha band [start_freq, end_freq]
%   beta_b      - Frequency range for the beta band [start_freq, end_freq]
%   spindle_b   - Frequency range for the spindle band [start_freq, end_freq]
%
% Outputs:
%   absolute_power - A structure containing the absolute power for each frequency band
%   relative_power - A structure containing the relative power for each frequency band (as a fraction of total power)
%   freq           - The frequency values (0.5 to 40 Hz) for the power spectrum
%   psd            - Power spectral density values corresponding to the frequency values
%
% Example:
%   [absolute_power, relative_power, freq, psd] = calculate_power_spectrum(EEG_clean, 1, 1:10000, 500, 0.5, [1 4], [4 8], [8 12], [12 30], [12 15]);
%
% The function calculates the absolute and relative power for each frequency band
% using the trapezoidal numerical integration (trapz) on the power spectral density (PSD).

function [absolute_power, relative_power, freq, psd] = calculate_power_spectrum(EEG_clean, chan, indexs, srate, freq_res, delta_b, theta_b, alpha_b, beta_b, spindle_b)
    % Function to calculate absolute and relative power in different frequency bands using pwelch.
    % Additionally, it outputs the frequency (freq) and power spectral density (psd).

    % Set up parameters for pwelch
    duration_window = 4 * srate;           % 4-second windows
    which_window = hanning(duration_window); % Hanning window
    overlap_window = duration_window * 0.5;    % 50% overlap (2 seconds)

    % Compute the power spectral density (psd) using pwelch
    [psd, freq] = pwelch(EEG_clean.data(chan, indexs), which_window, overlap_window, srate/freq_res, srate);

    % Filter frequencies of interest (0.5 to 40 Hz)
    freq_of_interest = freq >= 0.5 & freq <= 40;
    psd_filtered = psd(freq_of_interest, :);
    freq_filtered = freq(freq_of_interest, :);

    % Calculate absolute and relative power for each frequency band
    bands = {delta_b, theta_b, alpha_b, beta_b, spindle_b};
    band_names = {'delta', 'theta', 'alpha', 'beta', 'spindle'};
    absolute_power = struct();
    relative_power = struct();

    % Loop through each frequency band to calculate absolute and relative power
    for i = 1:length(bands)
        band = bands{i};
        band_name = band_names{i};

        % Indices corresponding to the current frequency band
        band_idx = freq_filtered >= band(1) & freq_filtered <= band(2);

        % Absolute power in the band
        absolute_power.(band_name) = trapz(freq_filtered(band_idx), psd_filtered(band_idx));

        % Relative power in the band (as a fraction of total power)
        relative_power.(band_name) = trapz(freq_filtered(band_idx), psd_filtered(band_idx)) / trapz(freq_filtered, psd_filtered);
    end

    % Return freq and psd for further analysis
    freq = freq_filtered;
    psd = psd_filtered;
end  
