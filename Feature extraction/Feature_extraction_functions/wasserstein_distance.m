% Function: wasserstein_distance
%
% Description:
% This function computes the Wasserstein distance between two EEG channels based on their power spectral density (PSD). 
% The PSD is calculated using the Welch method over the specified window and overlap parameters. The function then calculates 
% the cumulative distribution of the PSD in the frequency range of 0.5 to 40 Hz and compares the two channels using the 
% Wasserstein distance metric.
%
% If either of the channels is marked for removal or does not exist in the dataset, the function returns NaN for the 
% Wasserstein distance. The result is stored in a variable named dynamically based on the input channel names (e.g., 
% 'wasserstein_p3p4' for channels 'P3' and 'P4') and is displayed in the workspace.
%
% Inputs:
% - EEG_clean: Struct containing EEG data and channel information.
% - chan1: Name of the first channel (e.g., 'P3').
% - chan2: Name of the second channel (e.g., 'P4').
% - remove_channel_which: A logical array indicating channels to be removed from analysis.
% - which_window: The window function used for the Welch method.
% - overlap_window: The number of samples of overlap between adjacent segments.
% - duration_window: The length of the FFT window in samples.
% - srate: The sampling rate of the EEG data in Hz.
%
% Outputs:
% - The function dynamically creates a variable in the workspace named after the input channels (e.g., 'wasserstein_p3p4') 
%   and assigns it the Wasserstein distance value. The result is also displayed.
%
% Usage Example:
% wasserstein_distance(EEG_clean, 'P3', 'P4', remove_channel_which, hamming(128), 64, 256, 128);
%
% This function assumes that the input EEG dataset contains the specified channels and that proper window parameters are provided.




function wasserstein_distance(EEG_clean, chan1, chan2, remove_channel_which, which_window, overlap_window, duration_window, srate)
    % Check if either of the specified channels is marked for removal
    if remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan1)) || remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan2))
        wasserstein_dist = NaN;
    else
        % Ensure both channels exist in the dataset
        if all(~strcmp({EEG_clean.chanlocs.labels}, chan1)) || all(~strcmp({EEG_clean.chanlocs.labels}, chan2))
            wasserstein_dist = NaN;
        else
            % Compute the power spectral density (PSD) for both channels using pwelch
            [psd1, freq1] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan1), 1:length(EEG_clean.data)), ...
                which_window, overlap_window, duration_window, srate);
            [psd2, freq2] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan2), 1:length(EEG_clean.data)), ...
                which_window, overlap_window, duration_window, srate);

            % Calculate the cumulative distribution of the power spectral density in the 0.5 to 40 Hz range
            cum_psd1 = cumtrapz(psd1(freq1 >= 0.5 & freq1 <= 40));
            cum_psd1 = cum_psd1 ./ cum_psd1(end);
            cum_psd2 = cumtrapz(psd2(freq2 >= 0.5 & freq2 <= 40));
            cum_psd2 = cum_psd2 ./ cum_psd2(end);

            % Compute the Wasserstein distance between the two channels
            wasserstein_dist = sum(abs(cum_psd1 - cum_psd2)) / length(cum_psd1);
        end
    end
    
    % Create a variable name dynamically based on the input channel names
    variable_name = ['wasserstein_' lower(chan1) lower(chan2)];
    
    % Assign the Wasserstein distance to the workspace using the dynamically generated variable name
    assignin('base', variable_name, wasserstein_dist);
    
    % Display the result
    disp([variable_name, ' = ', num2str(wasserstein_dist)]);
end
