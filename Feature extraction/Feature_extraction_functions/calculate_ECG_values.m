function ECG_values = calculate_ECG_values(EEG_clean, indexs)
    % Function to analyze ECG data from EEG_clean and extract various HRV and complexity measures.
    % Inputs:
    %   - EEG_clean: EEG structure with ECG data in EEG_clean.ECG_data
    %   - indexs: Indices to extract the ECG data
    % Outputs:
    %   - ECG_values: A struct containing all calculated ECG measures
    
    % Sampling rate
    srate = EEG_clean.srate;
    
    % Initialize ECG_values struct
    ECG_values = struct();
    
    % Run Pan-Tompkins algorithm to detect QRS complexes
    [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(EEG_clean.data_ECG(1, indexs)), srate, 0);
    RR = diff(qrs_i_raw) / srate; % RR intervals in seconds
    
    if ~isempty(RR)
        ECG_values.ecg_mean = mean(60 ./ RR); % Heart rate in bpm
        
        % Normalize RR intervals
        RR_n = RR / mean(RR);
        ECG_values.ecg_quant20 = quantile(60 ./ RR_n, 0.2);
        ECG_values.ecg_quant80 = quantile(60 ./ RR_n, 0.8);
        
        % Time-domain measures
        ECG_values.ecg_rmssd = sqrt(mean(diff(RR_n) .^ 2));
        ECG_values.ecg_sdnn = std(RR_n);
        ECG_values.ecg_sdsd = std(diff(RR_n));
        ECG_values.ecg_pnn50 = sum(abs(diff(RR_n)) > 0.05) / length(RR) * 100; % Percentage of successive RR intervals differing by more than 50 ms
        
        % Frequency-domain measures
        Fs = 4; % Resampling frequency for frequency analysis
        t = cumsum(RR); % Time vector for RR intervals
        ti = 0:1/Fs:t(end); % Interpolated time vector
        RRi = interp1(t, RR, ti, 'pchip'); % Interpolate RR intervals
        
        % Power spectral density (PSD) using periodogram
        [Pxx, F] = periodogram(RRi, [], [], Fs);
        
        % Frequency bands
        ULF_band = F >= 0 & F <= 0.03;
        VLF_band = F > 0.03 & F <= 0.05;
        LF_band = F > 0.05 & F <= 0.15;
        HF_band = F > 0.15 & F <= 0.4;
        
        % Power in each band
        ECG_values.ecg_ULF = trapz(F(ULF_band), Pxx(ULF_band));
        ECG_values.ecg_VLF = trapz(F(VLF_band), Pxx(VLF_band));
        ECG_values.ecg_LF = trapz(F(LF_band), Pxx(LF_band));
        ECG_values.ecg_HF = trapz(F(HF_band), Pxx(HF_band));
        
        % Complexity measures (if there are enough RR intervals)
        if length(RR) > 10
            ECG_values.ecg_higuchi_fd = Higuchi_FD(zscore(RR), min(10, max(1, length(RR) - 2)));
            ECG_values.ecg_katz_fd = Katz_FD(zscore(RR));
            ECG_values.ecg_bubble_en = BubbEn(zscore(RR));
            ECG_values.ecg_spect_en = SpecEn(zscore(RR));
            ECG_values.ecg_sydy_en = SyDyEn(zscore(RR));
            ECG_values.ecg_phase_en = PhasEn(zscore(RR));
            ECG_values.ecg_slope_en = SlopEn(zscore(RR));
            ECG_values.ecg_eoe_en = EnofEn(zscore(RR));
            ECG_values.ecg_att_en = AttnEn(zscore(RR));
        else
            ECG_values.ecg_higuchi_fd = NaN;
            ECG_values.ecg_katz_fd = NaN;
            ECG_values.ecg_bubble_en = NaN;
            ECG_values.ecg_spect_en = NaN;
            ECG_values.ecg_sydy_en = NaN;
            ECG_values.ecg_phase_en = NaN;
            ECG_values.ecg_slope_en = NaN;
            ECG_values.ecg_eoe_en = NaN;
            ECG_values.ecg_att_en = NaN;
        end
    else
        % Set all ECG measures to NaN if no RR intervals are detected
        ECG_values.ecg_mean = NaN;
        ECG_values.ecg_quant20 = NaN;
        ECG_values.ecg_quant80 = NaN;
        ECG_values.ecg_rmssd = NaN;
        ECG_values.ecg_sdnn = NaN;
        ECG_values.ecg_sdsd = NaN;
        ECG_values.ecg_pnn50 = NaN;
        ECG_values.ecg_ULF = NaN;
        ECG_values.ecg_VLF = NaN;
        ECG_values.ecg_LF = NaN;
        ECG_values.ecg_HF = NaN;
        ECG_values.ecg_higuchi_fd = NaN;
        ECG_values.ecg_katz_fd = NaN;
        ECG_values.ecg_bubble_en = NaN;
        ECG_values.ecg_spect_en = NaN;
        ECG_values.ecg_sydy_en = NaN;
        ECG_values.ecg_phase_en = NaN;
        ECG_values.ecg_slope_en = NaN;
        ECG_values.ecg_eoe_en = NaN;
        ECG_values.ecg_att_en = NaN;
    end

    % ECG derived respiration rate (EDR)
    try
        py_peaks_result = py.neurokit2.ecg_peaks(py.numpy.array(EEG_clean.data_ECG(1, indexs)), ...
            sampling_rate=py.int(srate));
        py_ecg_rate_result = py.neurokit2.ecg_rate(py_peaks_result, sampling_rate=py.int(srate), desired_length=py.int(length(EEG_clean.data_ECG(1, indexs))));
        edr_result = double(py.neurokit2.ecg_rsp(py_ecg_rate_result, sampling_rate=srate).tolist());
        [edr_peak_amp, edr_peak_ind] = findpeaks(edr_result);
        ECG_values.br_permin = length(edr_peak_ind) / (length(indexs) / srate / 60);
    catch
        ECG_values.br_permin = NaN;
    end
end
