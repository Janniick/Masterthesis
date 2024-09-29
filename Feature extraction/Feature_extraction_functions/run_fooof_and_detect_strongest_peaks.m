function [fooof_offset, fooof_exponent, peak_table, normalized_absolute_power] = run_fooof_and_detect_strongest_peaks(freq, psd, delta_b, theta_b, alpha_b, beta_b, spindle_b, gamma_b)
    % Function to run FOOOF, detect the strongest peak in each frequency band, 
    % and calculate the normalized, absolute power for each band

    % Initialize FOOOF model
    f_result = py.fooof.FOOOF( ...
        peak_width_limits=[0.5, 12], ...    % Set peak width limits in Hz
        max_n_peaks=inf, ...                % Allow an unlimited number of peaks
        min_peak_height=0.0, ...            % Minimum peak height to consider
        peak_threshold=2.0, ...             % Threshold for peak detection
        aperiodic_mode='fixed', ...
        verbose=false);

    % Run FOOOF on the frequency spectrum
    f_result.fit(py.numpy.array(freq(freq >= delta_b(1) & freq <= gamma_b(2))'), ...
                 py.numpy.array(psd(freq >= delta_b(1) & freq <= gamma_b(2))'), ...
                 py.list([delta_b(1), gamma_b(2)]));

    % Extract aperiodic parameters (offset and exponent)
    aperpar = double(f_result.aperiodic_params_.tolist());
    fooof_offset = aperpar(1);  % Offset
    fooof_exponent = aperpar(2); % Exponent (1/f slope)

    % Extract cleaned spectrum (after removing 1/f)
    fooofed_spectrum = double(f_result.fooofed_spectrum_.tolist());
    ap_fit = py.fooof.sim.gen.gen_aperiodic(f_result.freqs, f_result.aperiodic_params_);
    normalized_spect = fooofed_spectrum - double(ap_fit.tolist());

    % Extract peak parameters
    outcell = cellfun(@double, cell(f_result.peak_params_.tolist()), 'UniformOutput', false);
    peak_params = reshape([outcell{:}], 3, size(outcell, 2))';

    % Define frequency bands and names
    bands = {delta_b, theta_b, alpha_b, beta_b, spindle_b, gamma_b};
    band_names = {'delta', 'theta', 'alpha', 'beta', 'spindle', 'gamma'};
    
    % Initialize the peak_table with a header row
    peak_table = {'Band', 'Frequency (Hz)', 'Amplitude', 'Bandwidth (Hz)'};
    normalized_absolute_power = struct();

    % Loop through bands to detect the strongest peak
    for i = 1:length(bands)
        band = bands{i};
        band_name = band_names{i};

        % Find peaks in the current frequency band
        band_idx = peak_params(:,1) >= band(1) & peak_params(:,1) < band(2);
        
        if any(band_idx)
            % Get the peaks in the current band
            peaks_in_band = peak_params(band_idx, :);
            
            % Find the strongest peak (based on amplitude)
            [~, strongest_idx] = max(peaks_in_band(:, 2));  % Find the peak with the maximum amplitude
            strongest_peak = peaks_in_band(strongest_idx, :);
            
            % Add the strongest peak to the table
            peak_table(end+1, :) = {band_name, strongest_peak(1), strongest_peak(2), strongest_peak(3)};
        else
            % No peak detected in this band
            peak_table(end+1, :) = {band_name, 'No peak', 'No peak', 'No peak'};
        end

        % Calculate normalized, absolute power for the current band
        normalized_absolute_power.(band_name) = sum(normalized_spect(freq >= band(1) & freq < band(2)));
    end

    % Calculate total power (for reference)
    normalized_absolute_power.total = sum(normalized_spect);

end
