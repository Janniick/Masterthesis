% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data';  % Replace with your folder path

% Search for all files ending with '_preprocessed.mat' in all subdirectories
all_mat_files = dir(fullfile(base_folder, '**', '*_preprocessed.mat'));

load([all_mat_files(i).folder, '\', all_mat_files(i).name]);

  %%%%% initialize things for spectrum and so on
% stfft sepctrogram
    which_EEG = EEG_clean;                   % which EEG signal
    srate = which_EEG.srate;                 % signal rate
    duration_window = 4*srate;               % 4 second windows
    which_window = hanning(duration_window); % hanning
    overlap_window = duration_window*0.5;    % 50%

    % multitaper spectrogram
    frequency_range = [0 40]; % Limit frequencies from 0 to 25 Hz
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

    %% Enlarge the scoring if it's condensed
    % Get the existing scoring table from EEG_clean
    scoring = EEG_clean.scoring;  % Assuming EEG_clean.scoring is the table
    % Initialize an empty table to store the expanded epochs
    expanded_scoring = table([], [], [], [], 'VariableNames', scoring.Properties.VariableNames);
    % Loop through each row of the original scoring table
    for i = 1:height(scoring)
        current_epoch = scoring.epoch(i);
        if i < height(scoring)
            next_epoch = scoring.epoch(i + 1);
        else
            next_epoch = current_epoch;  % Handle the last row
        end
        % Expand all epochs between current and next_epoch
        for epoch_num = current_epoch:(next_epoch - 1)
            new_row = table(epoch_num, scoring.stage(i), scoring.duration_s(i), scoring.sample_rate_hz(i), ...
                'VariableNames', scoring.Properties.VariableNames);
            expanded_scoring = [expanded_scoring; new_row];  % Append to the expanded table
        end
    end
    % If necessary, assign the expanded table back to EEG_clean.scoring
    EEG_clean.scoring_long = expanded_scoring;

    %% calculate how many fully clean epochs there are
    % Number of samples per epoch
    samples_per_epoch = EEG_clean.scoring_long.duration_s(1) * EEG_clean.scoring_long.sample_rate_hz(1);
    % Number of epochs in EEG_clean.scoring_long
    num_epochs = height(EEG_clean.scoring_long);
    % Initialize a counter for clean epochs
    clean_epochs_count = 0;
    % Loop through each epoch and check if it contains any bad signal (1's)
    for i = 1:num_epochs
        % Get the start and end indices for this epoch in EEG_clean.rem_ind_all
        start_idx = (i - 1) * samples_per_epoch + 1;
        end_idx = i * samples_per_epoch;
        % Extract the data for this epoch from rem_ind_all
        epoch_data = EEG_clean.rem_ind_all(start_idx:end_idx);
        % If there are no 1's in this epoch, count it as clean
        if all(epoch_data == 0)
            clean_epochs_count = clean_epochs_count + 1;
        end
    end
    % Calculate the percentage of clean epochs
    clean_epochs_percent = (clean_epochs_count / num_epochs) * 100;
    % Display the number of clean epochs
    disp(['Number of clean epochs: ', num2str(clean_epochs_count)]);
    disp(['Percentage of clean epochs: ', num2str(clean_epochs_percent), '%']);
    
    %% loop over channels
    nr_chan = length(which_EEG.chanlocs);
    %spect_channel = cell(nr_chan, 1);  % initialize channel data to save into full data
    eye_event_sample = eye_event*which_EEG.srate;
    if nr_chan > 0
        for chan = 1:nr_chan
            % fulldata = which_EEG.data(chan,:); % data to go through
            duration_fft = duration_window;    % duration_window or higher length(fulldata)
         
            %%%%%%%%%%%% cut algorithm
            new_cuts_idx = eeg_cut_and_glue(which_EEG, cut_chunks, chan, 0.8);


            %%%%%%%%%%%% main loops
            for eyes = ["open"] %["closed", "open"]
                if eyes == "closed"
                    start_index_o = 30*which_EEG.srate;
                    end_index_o = round(eye_event_sample-1);
                else
                    start_index_o = round(eye_event_sample) + 30*which_EEG.srate;
                    end_index_o = length(which_EEG.data(chan,:));
                end


                %%%%%%%%%% cutting out
                % cut out movements
                if ~isempty(new_cuts_idx)
                    % start to end, but cut out
                    index_cut = ismember(1:length(which_EEG.data(chan,:)), start_index_o:end_index_o);
                    for gii = 1:size(new_cuts_idx,1)
                        index_cut(new_cuts_idx(gii,1):new_cuts_idx(gii,2)) = false;
                    end
                else
                    index_cut = ismember(1:length(which_EEG.data(chan,:)), start_index_o:end_index_o);
                end

                % index on new scale
                index_cumsum = cumsum(index_cut);

                if sum(index_cut) == 0 
                    % no indices
                else
                    % calculate nr of epochs
                    full_len_s = max((length(which_EEG.data(chan,index_cut))/which_EEG.srate-30),0); % (length(which_EEG.data(chan,start_index_o:end_index_o))/which_EEG.srate-30);
                    nr_epochs = max(floor((full_len_s/epoch_windowlen - epoch_windowover)/(1-epoch_windowover)),0); % max(floor(full_len_s/epoch_windowlen),1);
    
                    if nr_epochs > 0
                        for epochs = 1:nr_epochs
                            start_index = (epochs-1)*epoch_windowlen*(1-epoch_windowover)*which_EEG.srate + 1; %(epochs-1)*epoch_windowlen*(1-epoch_windowover)*which_EEG.srate + start_index_o; % (epochs-1)*epoch_windowlen*which_EEG.srate + start_index_o;
                            end_index = min(start_index + epoch_windowlen*which_EEG.srate, max(index_cumsum)); % min(start_index + epoch_windowlen*which_EEG.srate, end_index_o) % min(epochs*epoch_windowlen*which_EEG.srate + start_index_o, end_index_o);
        
                            % previously everything from start_index:end_index
                            indexs = find(ismember(index_cumsum,start_index:end_index) & index_cut); % start_index:end:index
        
                            % bands #new
                            delta_ts = bandpass(which_EEG.data(chan, indexs), delta_b, which_EEG.srate);
                            theta_ts = bandpass(which_EEG.data(chan, indexs), theta_b, which_EEG.srate);
                            alpha_ts = bandpass(which_EEG.data(chan, indexs), alpha_b, which_EEG.srate);
                            beta_ts = bandpass(which_EEG.data(chan, indexs), beta_b, which_EEG.srate);
                            gamma_ts = bandpass(which_EEG.data(chan, indexs), gamma_b, which_EEG.srate);


                            % burst parameters #new
                            t_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(which_EEG.data(chan, indexs)), ...
                                fs=which_EEG.srate, dual_thresh=[1, 1.5], f_range=theta_b);
                            d_t_b = diff([0, double(t_burst.tolist()), 0]);
                            if sum(abs(d_t_b))/2 > 2
                                midpp = (find(d_t_b == 1) + find(d_t_b == -1) - 1) / 2;
                                t_bursts_ssd = std(diff(midpp)/which_EEG.srate);
                            else
                                t_bursts_ssd = NaN;
                            end
                            t_burst_stat = py.neurodsp.burst.compute_burst_stats(t_burst, fs=which_EEG.srate);
                            t_burst_stat = struct(t_burst_stat);
                            t_burst_stat.n_bursts = double(t_burst_stat.n_bursts);

                            a_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(which_EEG.data(chan, indexs)), ...
                                fs=which_EEG.srate, dual_thresh=[1, 1.5], f_range=alpha_b);
                            d_a_b = diff([0, double(a_burst.tolist()), 0]);
                            if sum(abs(d_a_b))/2 > 2
                                midpp = (find(d_a_b == 1) + find(d_a_b == -1) - 1) / 2;
                                a_bursts_ssd = std(diff(midpp)/which_EEG.srate);
                            else
                                a_bursts_ssd = NaN;
                            end
                            a_burst_stat = py.neurodsp.burst.compute_burst_stats(a_burst, fs=which_EEG.srate);
                            a_burst_stat = struct(a_burst_stat);
                            a_burst_stat.n_bursts = double(a_burst_stat.n_bursts);

                            b_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(which_EEG.data(chan, indexs)), ...
                                fs=which_EEG.srate, dual_thresh=[0.9, 1.5], f_range=beta_b);
                            d_b_b = diff([0, double(b_burst.tolist()), 0]);
                            if sum(abs(d_b_b))/2 > 2
                                midpp = (find(d_b_b == 1) + find(d_b_b == -1) - 1) / 2;
                                b_bursts_ssd = std(diff(midpp)/which_EEG.srate);
                            else
                                b_bursts_ssd = NaN;
                            end
                            b_burst_stat = py.neurodsp.burst.compute_burst_stats(b_burst, fs=which_EEG.srate);
                            b_burst_stat = struct(b_burst_stat);
                            b_burst_stat.n_bursts = double(b_burst_stat.n_bursts);

                            % coherence within #new
                            t_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(which_EEG.data(chan, indexs)), fs=which_EEG.srate, freqs=theta_b);
                            a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(which_EEG.data(chan, indexs)), fs=which_EEG.srate, freqs=alpha_b);
                            t_a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(which_EEG.data(chan, indexs)), fs=which_EEG.srate, freqs=[theta_b(1), alpha_b(2)]);
                            b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(which_EEG.data(chan, indexs)), fs=which_EEG.srate, freqs=beta_b);
                            a_b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(which_EEG.data(chan, indexs)), fs=which_EEG.srate, freqs=[alpha_b(1), beta_b(2)]);


           
                            % hurst long
                            hurs_exp_long = hurst_estimate(zscore(which_EEG.data(chan, indexs)), 'absval', 0, 1);
                            hurs_alpha = hurst_estimate(zscore(alpha_ts), 'absval', 0, 1);
                            hurs_beta = hurst_estimate(zscore(beta_ts), 'absval', 0, 1);
                            hurs_theta = hurst_estimate(zscore(theta_ts), 'absval', 0, 1);
                            hurs_delta = hurst_estimate(zscore(delta_ts), 'absval', 0, 1);
                             %#new
                            hurs_gamma = hurst_estimate(zscore(gamma_ts), 'absval', 0, 1);

                            % other fractal dimensions
                            higuchi_fd = Higuchi_FD(zscore(which_EEG.data(chan, indexs)), min(50, max(1,length(indexs)-2)));
                            katz_fd = Katz_FD(zscore(which_EEG.data(chan, indexs)));
                            bubble_en = BubbEn(zscore(which_EEG.data(chan, indexs)));
                            spect_en = SpecEn(zscore(which_EEG.data(chan, indexs)));
                            sydy_en = SyDyEn(zscore(which_EEG.data(chan, indexs)));
                            phase_en = PhasEn(zscore(which_EEG.data(chan, indexs)));
                            slope_en =  SlopEn(which_EEG.data(chan, indexs)); % original data
                            slope_en_z =  SlopEn(zscore(which_EEG.data(chan, indexs))); %zscore
                            slope_en_z_05_20 =  SlopEn(zscore(which_EEG.data(chan, indexs)), 'Lvls',[0.5,20]); %zscore
                            eoe_en = EnofEn(zscore(which_EEG.data(chan, indexs)));
                            att_en = AttnEn(zscore(which_EEG.data(chan, indexs)));
                             %#new
                            perm_en = PermEn(symb_king(which_EEG.data(chan, indexs)), 7);
                            perm_en = perm_en(7);
                            perm_en_t = PermEn(symb_king(theta_ts), 7);
                            perm_en_t = perm_en_t(7);
                            perm_en_a = PermEn(symb_king(alpha_ts), 7);
                            perm_en_a = perm_en_a(7);
                            perm_en_b = PermEn(symb_king(beta_ts), 7);
                            perm_en_b = perm_en_b(7);

                            % phase amp coupling #new
                            which_EEGt = which_EEG;
                            which_EEGt.data = which_EEGt.data(:,indexs);
                            evalc('which_EEGt = eeg_checkset(which_EEGt)');
                            which_EEGt = pop_pac(which_EEGt, 'channels', theta_b, alpha_b, chan,chan, 'forcecomp', 1); % pop_plotpac(EEG);
                            phase_amp_t_a = mean(which_EEGt.etc.eegpac(1).glm.pacval, 'all');
                            which_EEGt.etc.eegpac = [];
                            which_EEGt.etc.pacplotopt = [];
                            which_EEGt = pop_pac(which_EEGt, 'channels', alpha_b, beta_b, chan,chan, 'forcecomp', 1);
                            phase_amp_a_b = mean(which_EEGt.etc.eegpac(1).glm.pacval, 'all');
                            which_EEGt.etc.eegpac = [];
                            which_EEGt.etc.pacplotopt = [];
                            which_EEGt = pop_pac(which_EEGt, 'channels', theta_b, beta_b, chan,chan, 'forcecomp', 1);
                            phase_amp_t_b = mean(which_EEGt.etc.eegpac(1).glm.pacval, 'all');
                            which_EEGt.etc.eegpac = [];
                            which_EEGt.etc.pacplotopt = [];
                            clear which_EEGt;
            
                            % pwelch EOG
                            
                            if any(ismember(cell2table(use_neither_e12, 'VariableNames', {'c1', 'c2', 'c3'}), cell2table(current_infos, 'VariableNames', {'c1', 'c2', 'c3'}))) % test the one where both e1 and e2 are not usable
                                EOG.data = [];
                            end

                            if size(EOG.data,1)>1
                                duration_window = 2*srate;               % 2 second windows
                                which_window = hanning(duration_window); % hanning
                                overlap_window = duration_window*0.5;    % 50%
    
                                if any(ismember(cell2table(use_e1, 'VariableNames', {'c1', 'c2', 'c3'}), cell2table(current_infos, 'VariableNames', {'c1', 'c2', 'c3'}))) % test if current infos is in use_e1
                                    [psde1, freqe1] = pwelch(EOG.data(strcmp({EOG.chanlocs.labels},'E1'), indexs), ...
                                        which_window, overlap_window, duration_fft, srate);
                                    [psde2, freqe2] = pwelch(EEG.data(strcmp({EEG.chanlocs.labels},'C3'), indexs), ...
                                        which_window, overlap_window, duration_fft, srate);

                                    EOGd = EOG.data(strcmp({EOG.chanlocs.labels},'E1'),indexs);
                                    which_eog = 'E1';
                                else
                                    [psde1, freqe1] = pwelch(EOG.data(strcmp({EOG.chanlocs.labels},'E2'), indexs), ...
                                        which_window, overlap_window, duration_fft, srate);
                                    [psde2, freqe2] = pwelch(EEG.data(strcmp({EEG.chanlocs.labels},'C4'), indexs), ...
                                        which_window, overlap_window, duration_fft, srate);

                                    EOGd = EOG.data(strcmp({EOG.chanlocs.labels},'E2'),indexs);
                                    which_eog = 'E2';
                                end
                                
                                 %#new
                                id_all = find(freqe1 >= 0.5 & freqe1 <= 40);
                                id_swa = find(freqe1 >= delta_b(1) & freqe1 <= delta_b(2));
                                eog_d_abs = trapz(freqe1(id_swa), psde1(id_swa));
                                eog_d_rel = trapz(freqe1(id_swa), psde1(id_swa))/trapz(freqe1(id_all),psde1(id_all));

                                id_t = find(freqe1 >= theta_b(1) & freqe1 <= theta_b(2));
                                eog_ta_abs = trapz(freqe1(id_t), psde1(id_t));
                                eog_ta_rel = trapz(freqe1(id_t), psde1(id_t))/trapz(freqe1(id_all),psde1(id_all));

                                eog_c_d_abs = eog_d_abs - trapz(freqe2(id_swa), psde2(id_swa));
                                eog_c_d_rel = eog_d_rel - trapz(freqe1(id_swa), psde1(id_swa) - psde2(id_swa))/trapz(freqe2(id_all),psde1(id_all) - psde2(id_all));

                                id_ta = find(freqe1 >= theta_b(1) & freqe1 <= alpha_b(2));
                                eog_c_ta_abs = eog_ta_abs - trapz(freqe2(id_ta), psde2(id_ta));
                                eog_c_ta_rel = eog_ta_rel - trapz(freqe1(id_ta), psde1(id_ta) - psde2(id_ta))/trapz(freqe2(id_all),psde1(id_all) - psde2(id_all));

        
                                eog_higuchi_fd = Higuchi_FD(EOGd, min(50, max(1,length(indexs)-2)));
                                eog_katz_fd = Katz_FD(EOGd);
                                eog_hurst = hurst_estimate(EOGd, 'absval', 0, 1);
                                eog_slope = SlopEn(zscore(EOGd));

                                
                                
                                % Call the ecg_simulate function
                                py_eog_result = py.neurokit2.eog_process(EOGd, EOG.srate);
                                p2 = cell(py_eog_result);
                                p3 = struct(p2{2});
                                % p4 = table(p2{1});
                                eog_blinks_ind = length(int64(p3.EOG_Blinks.tolist()))+1;
                                % figure;plot(EOG.data(2,indexs)); hold on; xline(eog_blinks_ind);

                               



                            else
                                eog_d_abs = NaN;
                                eog_d_rel = NaN;
                                eog_ta_abs = NaN;
                                eog_ta_rel = NaN;

                                eog_c_d_abs = NaN;
                                eog_c_d_rel = NaN;
                                eog_c_ta_abs = NaN;
                                eog_c_ta_rel = NaN;

                                eog_blinks_ind = NaN;


                                eog_higuchi_fd = NaN;
                                eog_katz_fd = NaN;
                                eog_hurst = NaN;
                                eog_slope = NaN;
                            end
        
                    
                        
                        
            
                            % power density of current epoch
                            [psd, freq] = pwelch(which_EEG.data(chan, indexs), which_window, overlap_window, duration_fft, srate);

                            %#newnew timeseries alpha theta
                            [Swr, Fwr, ~] = spectrogram(which_EEG.data(chan, indexs), 4*srate, 2*srate, duration_fft, srate);
                            alpha_idx = find(Fwr >= alpha_b(1) & Fwr <= alpha_b(2));
                            theta_idx = find(Fwr >= theta_b(1) & Fwr <= theta_b(2));
                            alpha_power = sum(abs(Swr(alpha_idx, :)).^2, 1)*(srate/duration_fft);
                            theta_power = sum(abs(Swr(theta_idx, :)).^2, 1)*(srate/duration_fft);
                            q_ta = quantile(theta_power./alpha_power, [0.2, 0.8, 0.9, 1]);
                            q_ta_diff2080_20 = (q_ta(2)-q_ta(1))*q_ta(1);
                            q_ta_diff2080 = (q_ta(2)-q_ta(1));
                            q_ta_diff2090_20 = (q_ta(3)-q_ta(1))*q_ta(1);
                            q_ta_diff2090 = (q_ta(3)-q_ta(1));
                            q_ta_diff20max_20 = (q_ta(4)-q_ta(1))*q_ta(1);
                            q_ta_diff20max = (q_ta(4)-q_ta(1));
                            

                            
                    
                              %#new
                            % absolute and relative SWA of current epoch
                            id_swa = find(freq >= delta_b(1) & freq <= delta_b(2));
                            swa_abs = trapz(freq(id_swa), psd(id_swa));                       % similar as sum(psd(id_swa))*(srate/duration_fft)
                            swa_rel = trapz(freq(id_swa), psd(id_swa))/trapz(freq,psd);       % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));

                            % absolute and relative Theta activity
                            id_swa = find(freq >= theta_b(1) & freq <= theta_b(2));
                            ta_abs = trapz(freq(id_swa), psd(id_swa));
                            ta_rel = trapz(freq(id_swa), psd(id_swa))/trapz(freq,psd);        % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));

                            % absolute and relative Alpha activity
                            id_swa = find(freq >= alpha_b(1) & freq <= alpha_b(2));
                            aa_abs = trapz(freq(id_swa), psd(id_swa));
                            aa_rel = trapz(freq(id_swa), psd(id_swa))/trapz(freq,psd);        % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));

                            % absolute and relative Beta activity
                            id_swa = find(freq >= beta_b(1) & freq <= beta_b(2));
                            ba_abs = trapz(freq(id_swa), psd(id_swa));
                            ba_rel = trapz(freq(id_swa), psd(id_swa))/trapz(freq,psd);        % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));

                            %#newnew
                            % absolute and relative Gamma activity
                            id_swa = find(freq >= gamma_b(1) & freq <= gamma_b(2));
                            ga_abs = trapz(freq(id_swa), psd(id_swa));
                            ga_rel = trapz(freq(id_swa), psd(id_swa))/trapz(freq,psd);  
                            


                            % power density of occipital values
                            ochannels = ismember({which_EEG.chanlocs.labels}, {'O1', 'O2'});
                            if sum(ochannels)>0
                                eeg_ta_o_absolute = [];
                                eeg_ta_o_relative = [];
                                eeg_aa_o_absolute = [];
                                eeg_aa_o_relative = [];
                                for xo = find(ochannels)
                                    [psdo, freqo] = pwelch(which_EEG.data(xo, indexs), which_window, overlap_window, duration_fft, srate);

                                    %#new
                                    % absolute and relative Theta activity
                                    id_swa = find(freqo >= theta_b(1) & freqo <= theta_b(2));
                                    eeg_ta_o_absolute = [eeg_ta_o_absolute, trapz(freqo(id_swa), psdo(id_swa))];                          % similar as sum(psdo(id_swa))*(srate/duration_fft)
                                    eeg_ta_o_relative = [eeg_ta_o_relative, trapz(freqo(id_swa), psdo(id_swa))/trapz(freqo,psdo)];        % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));

                                    % absolute and relative Alpha activity
                                    id_swa = find(freqo >= alpha_b(1) & freqo <= alpha_b(2));
                                    eeg_aa_o_absolute = [eeg_aa_o_absolute, trapz(freqo(id_swa), psdo(id_swa))];
                                    eeg_aa_o_relative = [eeg_aa_o_relative, trapz(freqo(id_swa), psdo(id_swa))/trapz(freqo,psdo)];        % id_all = find(freq >= 0 & freq <= 40); trapz(freq(id_swa), psd(id_swa))/trapz(freq(id_all),psd(id_all));
                                end
                                eeg_ta_o_absolute = mean(eeg_ta_o_absolute);
                                eeg_ta_o_relative = mean(eeg_ta_o_relative);
                                eeg_aa_o_absolute = mean(eeg_aa_o_absolute);
                                eeg_aa_o_relative = mean(eeg_aa_o_relative);
                            else
                                eeg_ta_o_absolute = NaN;
                                eeg_ta_o_relative = NaN;
                                eeg_aa_o_absolute = NaN;
                                eeg_aa_o_relative = NaN;
                            end

                            

                            % other time vars
                            kurt = kurtosis(which_EEG.data(chan, indexs));
                            skew = skewness(which_EEG.data(chan, indexs));
                            [hjorth_activity, hjorth_mobility, hjorth_complexity] = hjorth(which_EEG.data(chan, indexs)',0);
            
                            % other complexity measures
            %                 which_EEG_down = pop_resample(which_EEG, 256);
            %                 [MSx,smse] = MSEn(which_EEG_down.data(chan, indexs), MSobject('SampEn')); % multiscale sample enntropy
            % 
            %                 scalz = unique(2*floor(linspace(which_EEG_down.srate/30, which_EEG_down.srate/0.7, 5)/2)+ 1);
            %                 mlz = getMultiscaleLZ(which_EEG_down.data(chan, indexs), scalz);
            %                 smlz = sum(mlz);
            
                            % mutual information in bands
                            mi_delta_theta = gcmi_cc(delta_ts(:), theta_ts(:));
                            mi_delta_alpha = gcmi_cc(delta_ts(:), alpha_ts(:));
                            mi_delta_beta = gcmi_cc(delta_ts(:), beta_ts(:));
                            mi_theta_alpha = gcmi_cc(theta_ts(:), alpha_ts(:));
                            mi_theta_beta = gcmi_cc(theta_ts(:), beta_ts(:));
                            mi_alpha_beta = gcmi_cc(alpha_ts(:), beta_ts(:));
            
                            % median mutual information between F and O channels

                                % anonymous if function hack
                                call = @(fun,par) fun(par{:});
                                cellget = @(cell, index) cell{index};
                                iff = @(test, truefn, falsefn, truePar, falsePar) call( ...
                                    cellget({falsefn, truefn}, test + 1), ... % functions
                                    cellget({falsePar, truePar}, test + 1) ... % params
                                );

                                checkzero = @(xy) iff(isempty(xy), @()false, @()xy == 0, {},{});

                                gcmi_cc_helper = @(chan_1, chan_2, remove_channel_which, which_EEG, indexs) ...
                                    iff(checkzero(remove_channel_which(strcmp({which_EEG.chanlocs.labels}, chan_1))) && ...
                                    checkzero(remove_channel_which(strcmp({which_EEG.chanlocs.labels}, chan_2))), ...
                                    @() gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels}, chan_1), indexs), ...
                                            which_EEG.data(strcmp({which_EEG.chanlocs.labels}, chan_2), indexs)), ...
                                    @() NaN, {}, {});

                                

                                mi_F_O_within = nanmean([gcmi_cc_helper('O1', 'F3', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('O2', 'F4', remove_channel_which, which_EEG, indexs)]);
                                mi_F_O_crossed = nanmean([gcmi_cc_helper('O1', 'F4', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('O2', 'F3', remove_channel_which, which_EEG, indexs)]);
                              
                                mi_F_C_within = nanmean([gcmi_cc_helper('C3', 'F3', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('C4', 'F4', remove_channel_which, which_EEG, indexs)]);
                                mi_F_C_crossed = nanmean([gcmi_cc_helper('C3', 'F4', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('C4', 'F3', remove_channel_which, which_EEG, indexs)]);

                                mi_C_O_within = nanmean([gcmi_cc_helper('O1', 'C3', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('O2', 'C4', remove_channel_which, which_EEG, indexs)]);
                                mi_C_O_crossed = nanmean([gcmi_cc_helper('O1', 'C4', remove_channel_which, which_EEG, indexs), ...
                                    gcmi_cc_helper('O2', 'C3', remove_channel_which, which_EEG, indexs)]);


                                mi_F_hemi = gcmi_cc_helper('F3', 'F4', remove_channel_which, which_EEG, indexs);
                                mi_C_hemi = gcmi_cc_helper('C3', 'C4', remove_channel_which, which_EEG, indexs);
                                mi_O_hemi = gcmi_cc_helper('O1', 'O2', remove_channel_which, which_EEG, indexs);

                                % try 
                                    % % very old version
                                % mi_F_O = median([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs))]);
                                % mi_F_C = median([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs))]);
                                % mi_C_O = median([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs))]);

                                % % old version
                                % mi_F_O_within = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs))]);
                                % mi_F_O_crossed = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs))]);
                                % 
                                % mi_F_C_within = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs))]);
                                % mi_F_C_crossed = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs))]);
                                % 
                                % mi_C_O_within = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs))]);
                                % mi_C_O_crossed = mean([gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs)), ...
                                %     gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs))]);
                                % 
                                % 
                                % mi_F_hemi = gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs));
                                % mi_C_hemi = gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs));
                                % mi_O_hemi = gcmi_cc(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                %     which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs));

                                % catch ME
                                %     mi_F_O_within = NaN;
                                %     mi_F_O_crossed = NaN;
                                % 
                                %     mi_F_C_within = NaN;
                                %     mi_F_C_crossed = NaN;
                                % 
                                %     mi_C_O_within = NaN;
                                %     mi_C_O_crossed = NaN;
                                % 
                                %     mi_F_hemi = NaN;
                                %     mi_C_hemi = NaN;
                                %     mi_O_hemi = NaN;
                                % end
            
                   
                            % fooof algorithm
                            f_result = py.fooof.FOOOF(peak_width_limits= [0.5 12], max_n_peaks=inf, min_peak_height=0.0, ...
                                peak_threshold=2.0, aperiodic_mode='fixed', verbose=false);
                            f_result.fit(py.numpy.array(freq'), py.numpy.array(psd'), py.list([1, 30]));

                            aperpar = double(f_result.aperiodic_params_.tolist());
                            fooofed_spectrum = double(f_result.fooofed_spectrum_.tolist());
                            ap_fit = py.fooof.sim.gen.gen_aperiodic(f_result.freqs, f_result.aperiodic_params_); %aperpar(1) - exp(aperpar(2))*foof_freq;
                            f_freqs = double(f_result.freqs.tolist());
                            f_exponent = aperpar(2);
                            f_offset = aperpar(1);
                            normalized_spect = fooofed_spectrum-double(ap_fit.tolist());
                            outcell = cellfun(@double,cell(f_result.peak_params_.tolist()),'UniformOutput',false);
                            peak_params = reshape([outcell{:}], 3, size(outcell,2))';

            
                            [mp, imp] = max(peak_params(:,2));
                            [mpa, impa] = max(peak_params(:,2).*peak_params(:,3));
                            
                            if isempty(imp)
                                strfpeak = NaN;
                                strfpeak_area = NaN;
                            else
                                strfpeak = peak_params(imp, 1);
                                strfpeak_area = peak_params(impa, 1);
                            end
            
                             %#new
                            d_sc = peak_params(:,1) >= delta_b(1) &  peak_params(:,1) < delta_b(2);
                            [mp, imp] = max(peak_params(d_sc,2));
                            [mpa, impa] = max(peak_params(d_sc,2).*peak_params(d_sc,3));
                            indf = find(d_sc);
                            if isempty(imp)
                                strfpeak_d = NaN;
                                strfpeak_d_area = NaN;
                            else
                                strfpeak_d = peak_params(indf(imp), 1);
                                strfpeak_d_area = peak_params(indf(impa), 1);
                            end

                            t_sc = peak_params(:,1) >= theta_b(1) &  peak_params(:,1) < theta_b(2);
                            [mp, imp] = max(peak_params(t_sc,2));
                            [mpa, impa] = max(peak_params(t_sc,2).*peak_params(t_sc,3));
                            indf = find(t_sc);
                            if isempty(imp)
                                strfpeak_t = NaN;
                                strfpeak_t_area = NaN;
                            else
                                strfpeak_t = peak_params(indf(imp), 1);
                                strfpeak_t_area = peak_params(indf(impa), 1);
                            end

                            t_sc = peak_params(:,1) >= theta_b(1)  &  peak_params(:,1) <= alpha_b(2);
                            [mp, imp] = max(peak_params(t_sc,2));
                            [mpa, impa] = max(peak_params(t_sc,2).*peak_params(t_sc,3));
                            indf = find(t_sc);
                            if isempty(imp)
                                strfpeak_ta = NaN;
                                strfpeak_ta_area = NaN;
                            else
                                strfpeak_ta = peak_params(indf(imp), 1);
                                strfpeak_ta_area = peak_params(indf(impa), 1);
                            end

                            a_sc = peak_params(:,1) >= alpha_b(1) &  peak_params(:,1) < alpha_b(2);
                            [mp, imp] = max(peak_params(a_sc,2));
                            [mpa, impa] = max(peak_params(a_sc,2).*peak_params(a_sc,3));
                            indf = find(a_sc);

                            if isempty(imp)
                                strfpeak_a = NaN;
                                strfpeak_a_area = NaN;
                            else
                                strfpeak_a = peak_params(indf(imp), 1);
                                strfpeak_a_area = peak_params(indf(impa), 1);
                            end

                            b_sc = peak_params(:,1) >= beta_b(1) &  peak_params(:,1) < beta_b(2);
                            [mp, imp] = max(peak_params(b_sc,2));
                            [mpa, impa] = max(peak_params(b_sc,2).*peak_params(b_sc,3));
                            indf = find(b_sc);
                            if isempty(imp)
                                strfpeak_b = NaN;
                                strfpeak_b_area = NaN;
                            else
                                strfpeak_b = peak_params(indf(imp), 1);
                                strfpeak_b_area = peak_params(indf(impa), 1);
                            end

                            %#new
                            g_sc = peak_params(:,1) >= gamma_b(1) &  peak_params(:,1) < gamma_b(2);
                            [mp, imp] = max(peak_params(g_sc,2));
                            [mpa, impa] = max(peak_params(g_sc,2).*peak_params(g_sc,3));
                            indf = find(g_sc);
                            if isempty(imp)
                                strfpeak_g = NaN;
                                strfpeak_g_area = NaN;
                            else
                                strfpeak_g = peak_params(indf(imp), 1);
                                strfpeak_g_area = peak_params(indf(impa), 1);
                            end


                            % periodic power #new
                            delta_f = sum(normalized_spect(f_freqs >= delta_b(1) & f_freqs < delta_b(2)));
                            theta_f = sum(normalized_spect(f_freqs >= theta_b(1) & f_freqs < theta_b(2)));
                            alpha_f = sum(normalized_spect(f_freqs >= alpha_b(1) & f_freqs < alpha_b(2)));
                            beta_f = sum(normalized_spect(f_freqs >= beta_b(1) & f_freqs < beta_b(2)));
                            %#new
                            gamma_f = sum(normalized_spect(f_freqs >= gamma_b(1) & f_freqs < gamma_b(2)));
                            total_f = sum(normalized_spect);
            
                            
                            % hemispheric spectral measure
                            
                            if remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'F4'))|remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'F3'))
                                wasserstein_frontal = NaN;
                            else
                                if all(~strcmp({which_EEG.chanlocs.labels}, 'F4'))|all(~strcmp({which_EEG.chanlocs.labels}, 'F3'))
                                    wasserstein_frontal = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F4'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'F3'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_frontal = sum(abs(cf3-cf4))/length(cf3);
                                end
                            end
            
                            if remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'C4'))|remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'C3'))
                                wasserstein_central = NaN;
                            else
                                if all(~strcmp({which_EEG.chanlocs.labels}, 'C4'))|all(~strcmp({which_EEG.chanlocs.labels}, 'C3'))
                                    wasserstein_central = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C4'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'C3'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_central = sum(abs(cf3-cf4))/length(cf3); 
                                end
                            end
            
                            if remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'O2'))|remove_channel_which(strcmp({which_EEG.chanlocs.labels}, 'O1')) 
                                wasserstein_occipital = NaN;
                            else
                                if all(~strcmp({which_EEG.chanlocs.labels}, 'O2'))|all(~strcmp({which_EEG.chanlocs.labels}, 'O1'))
                                    wasserstein_occipital = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O2'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(which_EEG.data(strcmp({which_EEG.chanlocs.labels},'O1'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_occipital = sum(abs(cf3-cf4))/length(cf3);
                                end
                            end
            
                            
            
            
            
                            % hilbert envelope hurst etc
                            hil = hilbert(which_EEG.data(chan, indexs));
                            hil_env = abs(hil);
            
                            hurst_hilbert = hurst_estimate(hil_env, 'absval', 0, 1);
                            slope_en_hilbert =  SlopEn(zscore(hil_env));
                            std_hilbert =  std(hil_env);
                            mean_hilbert =  mean(hil_env);
                            rmssd_hilbert = sqrt(mean(diff(hil_env).^2));
                            lzw_hilbert = lzw(hil_env>mean_hilbert);
            
        
        
        
                            % ECG measures
                            if isempty(ECG.data)|any(ismember(cell2table(remove_ecg, 'VariableNames', {'c1', 'c2', 'c3'}), cell2table(current_infos, 'VariableNames', {'c1', 'c2', 'c3'})))
                                ecg_mean = NaN;
                                ecg_quant20 = NaN;
                                ecg_quant80 = NaN;
                                ecg_rmssd = NaN;
                                ecg_sdnn = NaN;
                                ecg_sdsd = NaN;
                                ecg_pnn50 = NaN;
                                ecg_ULF = NaN;
                                ecg_VLF = NaN;
                                ecg_LF = NaN;
                                ecg_HF = NaN;

                                ecg_higuchi_fd = NaN;
                                ecg_katz_fd = NaN;
                                ecg_bubble_en = NaN;
                                ecg_spect_en = NaN;
                                ecg_sydy_en = NaN;
                                ecg_phase_en = NaN;
                                ecg_slope_en = NaN;
                                ecg_eoe_en = NaN;
                                ecg_att_en = NaN;
                                br_permin = NaN;
                            else
                                [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(ECG.data(1,indexs)),ECG.srate,0);
                                RR = diff(qrs_i_raw)/ECG.srate; % RR in seconds
            
                                 if ~isempty(RR)
                                    ecg_mean = mean(60./RR);
                                    %ecg_median = median(60./RR);

                                    % normalize
                                    RR_n = RR/mean(RR);
                                    ecg_quant20 = quantile(60./RR_n, 0.2);
                                    ecg_quant80 = quantile(60./RR_n, 0.8);
                                    
                                    % Time domain measures
                                    ecg_rmssd = sqrt(mean(diff(RR_n).^2));
                                    ecg_sdnn = std(RR_n);
                                    ecg_sdsd = std(diff(RR_n));

                                    ecg_pnn50 = sum(abs(diff(RR_n)) > 0.05) / length(RR) * 100; % take out?
                                    
                                    % Frequency domain measures
                                    Fs = 4; % Sampling frequency
                                    t = cumsum(RR); % Time vector
                                    ti = 0:1/Fs:t(end); % Interpolated time vector
                                    RRi = interp1(t, RR, ti, 'pchip');
                                    
                                    % Perform frequency analysis
                                    [Pxx, F] = periodogram(RRi, [], [], Fs);
                                    
                                    % Calculate power in each band
                                    % Define frequency bands
                                    ULF_band = F >= 0 & F <= 0.03;
                                    VLF_band = F > 0.03 & F <= 0.05;
                                    LF_band = F > 0.05 & F <= 0.15;
                                    HF_band = F > 0.15 & F <= 0.4;
                                    
                                    % Calculate power in each band using trapezoidal rule
                                    ecg_ULF = trapz(F(ULF_band), Pxx(ULF_band));
                                    ecg_VLF = trapz(F(VLF_band), Pxx(VLF_band));
                                    ecg_LF = trapz(F(LF_band), Pxx(LF_band));
                                    ecg_HF = trapz(F(HF_band), Pxx(HF_band));
                                    
                                    
                                    % Complexity measures
                                    if length(RR)>10
                                        ecg_higuchi_fd = Higuchi_FD(zscore(RR), min(10, max(1,length(RR)-2)));
                                        ecg_katz_fd = Katz_FD(zscore(RR));
                                        ecg_bubble_en = BubbEn(zscore(RR));
                                        ecg_spect_en = SpecEn(zscore(RR));
                                        ecg_sydy_en = SyDyEn(zscore(RR));
                                        ecg_phase_en = PhasEn(zscore(RR));
                                        ecg_slope_en =  SlopEn(zscore(RR));
                                        ecg_eoe_en = EnofEn(zscore(RR));
                                        ecg_att_en = AttnEn(zscore(RR));
                                    else
                                        ecg_higuchi_fd = NaN;
                                        ecg_katz_fd = NaN;
                                        ecg_bubble_en = NaN;
                                        ecg_spect_en = NaN;
                                        ecg_sydy_en = NaN;
                                        ecg_phase_en = NaN;
                                        ecg_slope_en =  NaN;
                                        ecg_eoe_en = NaN;
                                        ecg_att_en = NaN;
                                    end
                                else
                                    ecg_mean = NaN;
                                    ecg_quant20 = NaN;
                                    ecg_quant80 = NaN;
                                    ecg_rmssd = NaN;
                                    ecg_sdnn = NaN;
                                    ecg_sdsd = NaN;
                                    ecg_pnn50 = NaN;
                                    ecg_ULF = NaN;
                                    ecg_VLF = NaN;
                                    ecg_LF = NaN;
                                    ecg_HF = NaN;
                                    ecg_higuchi_fd = NaN;
                                    ecg_katz_fd = NaN;
                                    ecg_bubble_en = NaN;
                                    ecg_spect_en = NaN;
                                    ecg_sydy_en = NaN;
                                    ecg_phase_en = NaN;
                                    ecg_slope_en =  NaN;
                                    ecg_eoe_en = NaN;
                                    ecg_att_en = NaN;
                                end


                                % ecg derived respiration rate
                                try
                                    py_peaks_result = py.neurokit2.ecg_peaks(py.numpy.array(ECG.data(1,indexs)), ...
                                    sampling_rate=py.int(ECG.srate));
                                    py_ecg_rate_result = py.neurokit2.ecg_rate(py_peaks_result, sampling_rate=py.int(ECG.srate), desired_length=py.int(length(ECG.data(1,indexs))));
                                    edr_result = double(py.neurokit2.ecg_rsp(py_ecg_rate_result, sampling_rate=ECG.srate).tolist());
                                    [edr_peak_amp, edr_peak_ind] = findpeaks(edr_result);
                                    br_permin = length(edr_peak_ind)/(length(indexs)/ECG.srate/60);
                                catch ME
                                    br_permin = NaN;
                                end
                                
                            end
                           


                            

                            % breathing rate
                            % br_dat = ECG.data(1, min(indexs):max(indexs));
                            % if size(br_dat,1) < 1
                            %     br_mean = NaN;
                            %     br_rmssd = NaN;
                            %     br_sdnn = NaN;
                            %     br_sdsd = NaN;
                            % else
                            %     [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(br_dat),which_EEG.srate,0);
                            % 
                            %     % analysis_biosig = RunBioSigKit(br_dat, which_EEG.srate,0);
                            %     % edr = analysis_biosig.EDR_comp;
                            %     % sedr = smoothdata(edr, 'movmean',500);
                            %     % figure;plot(edr);hold on; plot(sedr);
                            %     % [s_pks,s_locs] = findpeaks(sedr);
                            %     % plot(s_locs, s_pks,'ob');
                            % 
                            %     rs_dist = zeros(1,length(qrs_i_raw));
                            %     s_peak = zeros(1,length(qrs_i_raw));
                            %     q_peak = zeros(1,length(qrs_i_raw));
                            %     for peaki = 1:length(qrs_i_raw)
                            %         pk_s = 0.05; % 50ms +-
                            %         window_pk = floor(pk_s*which_EEG.srate);
                            %         search_int = (qrs_i_raw(peaki) - window_pk):(qrs_i_raw(peaki) + window_pk);
                            %         look_pk = br_dat(search_int);
                            %         [mpk, ipk] = max(look_pk);
                            % 
                            %         % new peak aligned to actual maximum
                            %         qrs_i_raw(peaki) = search_int(ipk);
                            % 
                            %         % find S valley, calculate distance to
                            %         % R peak
                            %         search_int_s = (qrs_i_raw(peaki)+1):(qrs_i_raw(peaki)+window_pk);
                            %         look_s = br_dat(search_int_s);
                            %         [m_s, i_s] = min(look_s);
                            % 
                            %         s_peak(peaki) = search_int_s(i_s);
                            %         rs_dist(peaki) = br_dat(qrs_i_raw(peaki)) - br_dat(s_peak(peaki));
                            % 
                            %         % find Q valley
                            %         search_int_q = (qrs_i_raw(peaki)-window_pk):(qrs_i_raw(peaki)-1);
                            %         look_q = br_dat(search_int_q);
                            %         [m_q, i_q] = min(look_q);
                            % 
                            %         q_peak(peaki) = search_int_q(i_q);
                            %     end
                            % 
                            % 
                            %     s_fun = spline(qrs_i_raw, br_dat(qrs_i_raw));
                            %     s_x = min(qrs_i_raw):max(qrs_i_raw);
                            %     s_y = ppval(s_fun, s_x);
                            % 
                            % 
                            %     % lowpass filtering
                            %     mean_sy = mean(s_y);
                            %     padded_sy = [mean_sy*ones(1, which_EEG.srate*10), s_y, mean_sy*ones(1, which_EEG.srate*10)];
                            %     filtered_sy = lowpass(padded_sy, 0.5, which_EEG.srate);
                            %     filtered_sy = filtered_sy((which_EEG.srate*10+1):(end-which_EEG.srate*10));
                            % 
                            %     [s_pks,s_locs] = findpeaks(filtered_sy);
                            % 
                            %     figure;
                            %     plot(br_dat); hold on;
                            %     plot(qrs_i_raw, br_dat(qrs_i_raw)); 
                            %     plot(s_x, s_y)
                            %     plot(s_x, filtered_sy)
                            % 
                            %     BR = 60./(diff(s_locs)/which_EEG.srate);
                            % 
                            %     if length(BR) > 3
                            %         br_mean = mean(BR);
                            %         br_rmssd = sqrt(mean(diff(BR).^2));
                            %         br_sdnn = std(BR);
                            %         br_sdsd = std(diff(BR));
                            %     else
                            %         br_mean = NaN;
                            %         br_rmssd = NaN;
                            %         br_sdnn = NaN;
                            %         br_sdsd = NaN;
                            %     end
                            % end
              
            
            
        
                        
            
                            % save data #new
                            swa_data = [swa_data; table(string(thename), string(which_EEG.chanlocs(chan).labels), ...
                                eyes, epochs, ...
                                swa_abs, swa_rel, ...
                                ta_abs, ta_rel, ...
                                aa_abs, aa_rel, ...
                                ba_abs, ba_rel, ...
                                ga_abs, ga_rel, ...
                                eeg_ta_o_absolute, ...
                                eeg_ta_o_relative, ...
                                eeg_aa_o_absolute, ...
                                eeg_aa_o_relative, ...
                                wasserstein_frontal, wasserstein_central, wasserstein_occipital, ...
                                f_exponent, f_offset, ...
                                delta_f, theta_f, alpha_f, beta_f, gamma_f, total_f, ...
                                strfpeak, strfpeak_area, ...
                                strfpeak_d, strfpeak_d_area, ...
                                strfpeak_t, strfpeak_t_area, ...
                                strfpeak_ta, strfpeak_ta_area, ...
                                strfpeak_a, strfpeak_a_area, ...
                                strfpeak_b, strfpeak_b_area, ...
                                strfpeak_g, strfpeak_g_area, ...
                                eog_d_abs, ...
                                eog_d_rel, ...
                                eog_ta_abs, ...
                                eog_ta_rel, ...
                                eog_c_d_abs, ...
                                eog_c_d_rel, ...
                                eog_c_ta_abs, ...
                                eog_c_ta_rel, ...
                                eog_blinks_ind, ...
                                eog_higuchi_fd, ...
                                eog_katz_fd, ...
                                eog_hurst, ...
                                eog_slope, ...
                                hurs_exp_long, ...
                                hurs_alpha, ...
                                hurs_beta, ...
                                hurs_theta, ...
                                hurs_delta, ... %, blpm, bld, blpa, blna, ...
                                higuchi_fd, ...
                                katz_fd, ...
                                kurt, ...
                                skew, ...
                                hjorth_activity, ...
                                hjorth_mobility, ...
                                hjorth_complexity, ...
                                mi_delta_theta, ...
                                mi_delta_alpha, ...
                                mi_delta_beta, ...
                                mi_theta_alpha, ...
                                mi_theta_beta, ...
                                mi_alpha_beta, ...
                                mi_F_O_within, ...
                                mi_F_O_crossed, ...
                                mi_F_C_within, ...
                                mi_F_C_crossed, ...
                                mi_C_O_within, ...
                                mi_C_O_crossed, ...
                                mi_F_hemi, ...
                                mi_C_hemi, ...
                                mi_O_hemi, ...
                                bubble_en, ...
                                spect_en, ...
                                sydy_en, ...
                                phase_en, ...
                                slope_en, ...
                                slope_en_z, ...
                                slope_en_z_05_20, ...
                                eoe_en, ...
                                att_en, ...
                                hurst_hilbert, ...
                                slope_en_hilbert, ...
                                std_hilbert, ...
                                mean_hilbert, ...
                                rmssd_hilbert, ...
                                lzw_hilbert, ...
                                ecg_mean, ...
                                ecg_quant20, ...
                                ecg_quant80, ...
                                ecg_rmssd, ...
                                ecg_sdnn, ...
                                ecg_sdsd, ...
                                ecg_pnn50, ...
                                ecg_ULF, ...
                                ecg_VLF, ...
                                ecg_LF, ...
                                ecg_HF, ...
                                ecg_higuchi_fd, ...
                                ecg_katz_fd, ...
                                ecg_bubble_en, ...
                                ecg_spect_en, ...
                                ecg_sydy_en, ...
                                ecg_phase_en, ...
                                ecg_slope_en, ...
                                ecg_eoe_en, ...
                                ecg_att_en, ...
                                br_permin, ...
                                t_burst_stat.n_bursts, ...
                                t_burst_stat.duration_mean, ...
                                t_burst_stat.duration_std, ...
                                a_burst_stat.n_bursts, ...
                                a_burst_stat.duration_mean, ...
                                a_burst_stat.duration_std, ...
                                b_burst_stat.n_bursts, ...
                                b_burst_stat.duration_mean, ...
                                b_burst_stat.duration_std, ...
                                t_lagcoh, ...
                                a_lagcoh, ...
                                t_a_lagcoh, ...
                                b_lagcoh, ...
                                a_b_lagcoh, ...
                                t_bursts_ssd, ...
                                a_bursts_ssd, ...
                                b_bursts_ssd, ...
                                phase_amp_t_a, ...
                                phase_amp_a_b, ...
                                phase_amp_t_b, ...
                                perm_en, ...
                                perm_en_t, ...
                                perm_en_a, ...
                                perm_en_b, ...
                                q_ta_diff2080_20, ...
                                q_ta_diff2080, ...
                                q_ta_diff2090_20, ...
                                q_ta_diff2090, ...
                                q_ta_diff20max_20, ...
                                q_ta_diff20max, ...
                                'VariableNames', variable_names_types(:,1))];
                        end
                    end



%%%%% other script with parallel processing and powerspectrum and foof
 if size(EEG.data,1) > 0 & ~isempty(current_scoring)

        nr_epochs = min(ceil((size(EEG.data,2)/EEG.srate)/30),size(current_scoring,1));
        if ~ ceil((size(EEG.data,2)/EEG.srate)/30) == size(current_scoring,1)
            fprintf('------------Warning nr_epochs not matching for iteration %d\n', i);
        end
        nr_channels = size(EEG.data,1);


        % set up pwelch method for power spectrum
        duration_window = 4*EEG.srate;           % 4s windows
        which_window = hanning(duration_window); % hanning window
        overlap_window = duration_window*0.5;    % 50% overlap, 2s

        % setup epochs
        epoch_windowlen = 30; % seconds
        epoch_windowover = 0; % no overlap

        % if enough epochs
        if length(EEG.data(1, :))/EEG.srate >= 2 && nr_epochs >= 2
            % start cluster with 80% of cores
            numCores = feature('numcores');
            numWorkers = max(1, floor(0.8 * numCores));

            for chan = 1:nr_channels
                disp(EEG.chanlocs(chan).labels);
                temp_t = cell(nr_epochs, 1);

                % start cluster
                parpool(numWorkers);
                parfor epo = 1:nr_epochs
                    % fprintf('%d ',epo)
                    start_index = (epo-1)*epoch_windowlen*(1-epoch_windowover)*EEG.srate + 1;
                    end_index = min(start_index + epoch_windowlen*EEG.srate, size(EEG.data,2));
                    idxs = start_index:end_index;

                    per_badi = sum(rem_ind_all(idxs))/(epoch_windowlen*EEG.srate);
                    if per_badi <= 0.1 && (end_index - start_index)*(1-per_badi)> length(which_window)
                        idxs = idxs(~rem_ind_all(idxs)); % take only the good ones
                        % do pwelch spectrum
                        [psd, freq] = pwelch(EEG.data(chan,idxs), ...
                            which_window, overlap_window, EEG.srate/freq_res, EEG.srate);

                        freq_of_interest = freq >= 0.5 & freq <= 30;
                        psdl = 10*log10(psd(freq_of_interest,:))';
                        psdo = psd(freq_of_interest,:)';
                        freql = freq(freq_of_interest,:)';

                        id_swa = find(freql >= 0.5 & freql <= 2);
                        eeg_d1_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_d1_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 2 & freql <= 4);
                        eeg_d2_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_d2_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 0.5 & freql <= 4);
                        eeg_d_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_d_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 4 & freql <=6);
                        eeg_t1_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_t1_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 6 & freql <= 8);
                        eeg_t2_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_t2_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 4 & freql <= 8);
                        eeg_t_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_t_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 8 & freql <= 10);
                        eeg_a1_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_a1_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 10 & freql <= 12);
                        eeg_a2_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_a2_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 8 & freql <= 12);
                        eeg_a_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_a_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 13 & freql <= 20);
                        eeg_b1_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_b1_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >=20 & freql <= 30);
                        eeg_b2_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_b2_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 15 & freql <= 30);
                        eeg_b_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_b_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        id_swa = find(freql >= 11 & freql <= 13);
                        eeg_s1_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_s1_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 13 & freql <= 16);
                        eeg_s2_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_s2_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);
                        id_swa = find(freql >= 11 & freql <= 16);
                        eeg_s_abs = trapz(freql(id_swa), psdo(id_swa));
                        eeg_s_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                        % fooof algorithm
                        f_result = py.fooof.FOOOF(peak_width_limits= [0.5 12], max_n_peaks=inf, min_peak_height=0.0, ...
                            peak_threshold=2.0, aperiodic_mode='fixed', verbose=false);
                        f_result.fit(py.numpy.array(freq'), py.numpy.array(psd'), py.list([1, 30]));

                        aperpar = double(f_result.aperiodic_params_.tolist());
                        f_exponent = aperpar(2);
                        f_offset = aperpar(1);


                        % save
                        stagecurr = string(current_scoring.stage(current_scoring.epoch == epo));
                        if isempty(stagecurr)
                            stagecurr = NaN;
                        end
                        temp_t{epo} = cell2table([{string(thenamescore(2))}, ...
                            {string(thenamescore(3))}, ...
                            {epo}, ...
                            {stagecurr(1)}, ... % 'MGT_H23_N4_sleep_sienna2_808r0s2' i = 59 has two stages named 1, 2, 3,....
                            {chan}, ...
                            {f_exponent}, ...
                            {f_offset}, ...
                            num2cell([eeg_d1_abs, eeg_d1_rel, eeg_d2_abs, eeg_d2_rel, eeg_d_abs, eeg_d_rel,...
                            eeg_t1_abs, eeg_t1_rel, eeg_t2_abs, eeg_t2_rel, eeg_t_abs, eeg_t_rel,...
                            eeg_a1_abs, eeg_a1_rel, eeg_a2_abs, eeg_a2_rel, eeg_a_abs, eeg_a_rel,...
                            eeg_b1_abs, eeg_b1_rel, eeg_b2_abs, eeg_b2_rel, eeg_b_abs, eeg_b_rel,...
                            eeg_s1_abs, eeg_s1_rel, eeg_s2_abs, eeg_s2_rel, eeg_s_abs, eeg_s_rel]), ...
                            num2cell(psdo)], 'VariableNames', colnames_spect);
                    else % end if larger than window size
                        temp_t{epo} =[];
                    end
                end % epoch loop
                delete(gcp('nocreate'));

                temp_t = temp_t(~cellfun('isempty', temp_t));
                spectrum_info = [spectrum_info; vertcat(temp_t{:})];
            end % channel loop
        end % enough epochs
    end % enough data and scoring




###### again other script, with SW detection and travelling waves
if size(EEG.data,1) > 0 & ~isempty(current_scoring)
        nr_epochs = min(ceil((size(EEG.data,2)/EEG.srate)/30),size(current_scoring,1));
        if ~ ceil((size(EEG.data,2)/EEG.srate)/30) == size(current_scoring,1)
            fprintf('------------Warning nr_epochs not matching for iteration %d\n', i);
        end
        nr_channels = size(EEG.data,1);

     
        epoch_windowlen = 30; % seconds
        epoch_windowover = 0; % no overlap
        if length(EEG.data(1, :))/EEG.srate >= 2
            % add stage to each data sample
            current_stages = str2double(string(current_scoring.stage)');
            lastvalue = find(~isnan(table2array(current_scoring(:,'interval_s')))');
            lastvalue = lastvalue(end);
            hypnodata = [repelem(current_stages(1),table2array(current_scoring(1,'start_s'))*EEG.srate), ...
                        repelem(current_stages(1:lastvalue), 30*EEG.srate), ...
                        repelem(current_stages(lastvalue), max(0,(size(EEG.data,2)/EEG.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG.srate))
                        ];
            current_epoch = str2double(string(current_scoring.epoch)');
            epochdata = [repelem(current_epoch(1),table2array(current_scoring(1,'start_s'))*EEG.srate), ...
                        repelem(current_epoch(1:lastvalue), 30*EEG.srate), ...
                        repelem(current_epoch(lastvalue), max(0,(size(EEG.data,2)/EEG.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG.srate))
                        ];
                    
                                
            data_for_py = cell(1, size(EEG.data, 1));
            for row = 1:size(EEG.data, 1)
                data_for_py(row) = {EEG.data(row, :)};
            end
            data_for_py = py.numpy.array(data_for_py);


            sws_detection_result = py.yasa.sw_detect(data = data_for_py, ...
                        sf = py.float(EEG.srate), ...
                        ch_names=py.list({EEG.chanlocs.labels}), ...
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


            if ~isempty(sws_detection_info)
                sws_detection_info = sws_detection_result.summary();
                sws_colnames = cellfun(@string,cell(sws_detection_info.columns.values.tolist()));
                sws_detection_info = cell(sws_detection_info.values.tolist());
                sws_detection_info = cellfun(@(xs)cell2table(cell(xs), 'VariableNames', sws_colnames),sws_detection_info, 'UniformOutput', false);
                sws_detection_info = vertcat(sws_detection_info{:});
                sws_detection_info.Stage = cellfun(@double, sws_detection_info.Stage);
                sws_detection_info.Channel = string(sws_detection_info.Channel);
                sws_detection_info.IdxChannel = cellfun(@double, sws_detection_info.IdxChannel);

                sws_detection_info.neg_sw_peaks_ind = sws_detection_info.NegPeak*EEG.srate;
                sws_detection_info.duration_negative_phase = (sws_detection_info.MidCrossing - sws_detection_info.Start)*EEG.srate;
                sws_detection_info.duration_positive_phase = (sws_detection_info.End - sws_detection_info.MidCrossing)*EEG.srate;
                sws_detection_info.duration_start_to_neg = (sws_detection_info.NegPeak - sws_detection_info.Start)*EEG.srate;
                sws_detection_info.duration_neg_to_mid = (sws_detection_info.MidCrossing - sws_detection_info.NegPeak)*EEG.srate;
                sws_detection_info.duration_mid_to_pos = (sws_detection_info.PosPeak - sws_detection_info.MidCrossing)*EEG.srate;
                sws_detection_info.duration_pos_to_end = (sws_detection_info.End - sws_detection_info.PosPeak)*EEG.srate;


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
                        min_20ms = ref_neg_idx - floor(0.02*EEG.srate);
                        plus_20ms = ref_neg_idx + floor(0.02*EEG.srate);

                        sws_detection_info.percent_involved(idx_refwave(refw)) = sum(mean(EEG.data(:, min_20ms:plus_20ms),2)< -5)/size(EEG.data,1);

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


