function complexity_measures = extract_complexity_measures(EEG_clean, srate, delta_b, theta_b, alpha_b, beta_b, spindle_b)
    % extract_complexity_measures - Compute and organize EEG complexity measures.
    %
    % Syntax: complexity_measures = extract_complexity_measures(EEG_clean, srate, delta_b, theta_b, alpha_b, beta_b, spindle_b)
    %
    % Inputs:
    %    EEG_clean   - Struct containing EEG data and channel information.
    %    srate       - Sampling rate of EEG data in Hz.
    %    delta_b     - 1x2 array defining delta band [low_freq, high_freq].
    %    theta_b     - 1x2 array defining theta band [low_freq, high_freq].
    %    alpha_b     - 1x2 array defining alpha band [low_freq, high_freq].
    %    beta_b      - 1x2 array defining beta band [low_freq, high_freq].
    %    spindle_b   - 1x2 array defining spindle band [low_freq, high_freq].
    %
    % Outputs:
    %    complexity_measures - Struct containing all computed complexity measures.
    %
    % Author:
    %    Jannick
    %
    % Date:
    %    October 2024

    %% Initialize the Complexity Measures Struct
    complexity_measures = struct();
    
    %% Define Frequency Bands for Iteration
    bands = {'delta', 'theta', 'alpha', 'beta', 'spindle'};
    band_ranges = {delta_b, theta_b, alpha_b, beta_b, spindle_b};
    
    %% A. Full Night, Whole Data
    fprintf('Computing Full Night, Whole Data Complexity Measures...\n');
    complexity_measures.full_night.whole_data = struct();
    
    % 1. Hurst Exponent
    fprintf(' - Computing Hurst Exponent...\n');
    complexity_measures.full_night.whole_data.hurst = hurst_estimate(zscore(EEG_clean.data(:)), 'absval', 0, 1);
    
    % 2. Hjorth Parameters
    fprintf(' - Computing Hjorth Parameters...\n');
    [hjorth_activity, hjorth_mobility, hjorth_complexity] = hjorth(EEG_clean.data(:)', 0);
    complexity_measures.full_night.whole_data.hjorth_activity = hjorth_activity;
    complexity_measures.full_night.whole_data.hjorth_mobility = hjorth_mobility;
    complexity_measures.full_night.whole_data.hjorth_complexity = hjorth_complexity;
    
    % 3. Statistical Measures
    fprintf(' - Computing Statistical Measures (Kurtosis & Skewness)...\n');
    complexity_measures.full_night.whole_data.kurtosis = kurtosis(EEG_clean.data(:));
    complexity_measures.full_night.whole_data.skewness = skewness(EEG_clean.data(:));
    
    % 4. Slope Entropy
    fprintf(' - Computing Slope Entropy...\n');
    complexity_measures.full_night.whole_data.slope_entropy = SlopEn(zscore(EEG_clean.data(:)));
    
    % % 5. Multiscale Measures
    % dauert unendlich lange
    % fprintf(' - Computing Multiscale Sample Entropy and Multiscale Lempel-Ziv Complexity...\n');
    % %complexity_measures.full_night.whole_data.multiscale_sample_entropy = MSEn(EEG_clean.data(:), MSobject('SampEn'));
    % scalz = unique(2 * floor(linspace(srate/30, srate/0.7, 5)/2) + 1); 
    % complexity_measures.full_night.whole_data.multiscale_lz_complexity = sum(getMultiscaleLZ(EEG_clean.data(:), scalz));
    
    %% B. Full Night, Channel-Specific
    fprintf('Computing Full Night, Channel-Specific Complexity Measures...\n');
    num_channels = length(EEG_clean.chanlocs);
    channels = {EEG_clean.chanlocs.labels};
    
    % Initialize Channel-Specific Struct
    complexity_measures.full_night.channel_specific = struct();
    
    for ch = 1:num_channels
        channel_label = channels{ch};
        fprintf(' - Processing Channel: %s\n', channel_label);
        complexity_measures.full_night.channel_specific.(channel_label) = struct();
        
        % 1. Permutation Entropy
        % dauert sehr lange
        % fprintf('   * Computing Permutation Entropy...\n');
        % symb_signal = symb_king(EEG_clean.data(ch, :));
        % perm_en = PermEn(symb_signal, 7);
        % complexity_measures.full_night.channel_specific.(channel_label).permutation_entropy = perm_en(7);
        
        % 2. Spectral Entropy
        %dauert lange
        % fprintf('   * Computing Spectral Entropy...\n');
        % spect_en = SpecEn(zscore(EEG_clean.data(ch, :)));
        % complexity_measures.full_night.channel_specific.(channel_label).spectral_entropy = spect_en;
       
        % % 3. Fractal Dimensions
        % dauert lange
        % fprintf('   * Computing Fractal Dimensions (Katz)...\n');
        % complexity_measures.full_night.channel_specific.(channel_label).katz_fd = Katz_FD(zscore(EEG_clean.data(ch, :)));
        
        % 4. Hjorth Parameters
        fprintf('   * Computing Hjorth Parameters...\n');
        [hj_act, hj_mob, hj_comp] = hjorth(EEG_clean.data(ch, :)', 0);
        complexity_measures.full_night.channel_specific.(channel_label).hjorth_activity = hj_act;
        complexity_measures.full_night.channel_specific.(channel_label).hjorth_mobility = hj_mob;
        complexity_measures.full_night.channel_specific.(channel_label).hjorth_complexity = hj_comp;
        
        % 5. Slope Entropy with Specific Levels
        fprintf('   * Computing Slope Entropy with Levels [0.5, 20]...\n');
        complexity_measures.full_night.channel_specific.(channel_label).slope_entropy_z_05_20 = SlopEn(zscore(EEG_clean.data(ch, :)), 'Lvls', [0.5, 20]);
        
    end
    
    %% C. Full Night, Band-Specific
    fprintf('Computing Full Night, Band-Specific Complexity Measures...\n');
    complexity_measures.full_night.band_specific = struct();
    
    for b = 1:length(bands)
        band = bands{b};
        band_range = band_ranges{b};
        fprintf(' - Processing Frequency Band: %s [%d-%d] Hz\n', band, band_range(1), band_range(2));
        complexity_measures.full_night.band_specific.(band) = struct();
        
        % Extract the signal filtered for the current band
        fprintf('   * Filtering Signal for %s Band...\n', band);
        % Assuming you have a function or method to filter the EEG signal into bands
        band_signal = bandpass(EEG_clean.data, band_range, srate);
        
        % 1. Permutation Entropy per Band
        % dauert lange
        % fprintf('   * Computing Permutation Entropy for %s Band...\n', band);
        % perm_en_band = PermEn(symb_king(band_signal(:)), 7);
        % complexity_measures.full_night.band_specific.(band).permutation_entropy = perm_en_band(7);
        
        % 2. Spectral Entropy per Band
        % dauert lange
        % fprintf('   * Computing Spectral Entropy for %s Band...\n', band);
        % spect_en_band = SpecEn(zscore(band_signal(:)), 7);  % Assuming SpecEn accepts the entire band signal
        % complexity_measures.full_night.band_specific.(band).spectral_entropy = spect_en_band(7);
        
        % 3. Hurst Exponent per Band
        fprintf('   * Computing Hurst Exponent for %s Band...\n', band);
        complexity_measures.full_night.band_specific.(band).hurst_exponent = hurst_estimate(zscore(band_signal(:)), 'absval', 0, 1);
        
        % 4. Lagged Coherence per Band
        fprintf('   * Computing Lagged Coherence for %s Band...\n', band);
        % Convert band_signal to Python NumPy array
        band_signal_py = py.numpy.array(band_signal);
        % Compute lagged coherence using neurodsp.rhythm.compute_lagged_coherence
        try
            lag_coh = py.neurodsp.rhythm.compute_lagged_coherence(band_signal_py, pyargs('fs', srate, 'freqs', band_range));
            % Convert Python array to MATLAB array
            lag_coh_matlab = double(lag_coh);
            % Compute mean lagged coherence (if multiple values are returned)
            mean_lag_coh = mean(lag_coh_matlab);
            complexity_measures.full_night.band_specific.(band).lagged_coherence = mean_lag_coh;
        catch ME
            warning('Lagged Coherence computation failed for band %s: %s', band, ME.message);
            complexity_measures.full_night.band_specific.(band).lagged_coherence = NaN;
        end
    end
    
    %% Finalizing the Struct
    fprintf('Complexity Measures Computation Completed.\n');
    
    end
 
