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

            % Calculate the cumulative distribution of the power spectral density in the 0.5 to 30 Hz range
            cum_psd1 = cumtrapz(psd1(freq1 >= 0.5 & freq1 <= 30));
            cum_psd1 = cum_psd1 ./ cum_psd1(end);
            cum_psd2 = cumtrapz(psd2(freq2 >= 0.5 & freq2 <= 30));
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
