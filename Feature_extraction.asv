%% add functions
addpath(genpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Feature extraction\Feature_extraction_functions'));
addpath(genpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_mesmart\scripts'));
addpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\matlab_processing_benji');
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
%% load python
pyenv('Version', 'C:\Users\janni\AppData\Local\Programs\Python\Python311\python.EXE');
%% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data';  % Replace with your folder path

% Search for all files ending with '_preprocessed.mat' in all subdirectories
all_mat_files = dir(fullfile(base_folder, '**', '*_preprocessed.mat'));

load([all_mat_files(i).folder, '\', all_mat_files(i).name]);

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
    
    %% convert 20s to 30s scoring if needed
    EEG_clean.scoring_long_30s = scoring_convert20sTo30s(EEG_clean.scoring_long);
  %% Check for clean epochs
% Number of samples per epoch
samples_per_epoch = EEG_clean.scoring_long_30s.duration_s(1) * EEG_clean.scoring_long_30s.sample_rate_hz(1);
% Number of epochs in EEG_clean.scoring_long_30s
num_epochs = height(EEG_clean.scoring_long_30s);
% Initialize a counter for clean epochs and an array to store clean epoch data
clean_epochs_count = 0;
clean_epoch_data = [];  % Array to store [clean_epoch, start_idx, end_idx]
clean_epoch_info = [];  % To store the additional entries from scoring_long_30s

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
        % Store [clean_epoch, start_idx, end_idx]
        clean_epoch_data = [clean_epoch_data; i, start_idx, end_idx];
        % Also store the corresponding entries from scoring_long_30s
        clean_epoch_info = [clean_epoch_info; EEG_clean.scoring_long_30s(i, :)];
    end
end

% Convert the clean_epoch_data array to a table
clean_epochs_table = array2table(clean_epoch_data, 'VariableNames', {'clean_epoch', 'start_epoch_index', 'end_epoch_index'});

% Add the columns from EEG_clean.scoring_long_30s to the clean_epochs_table
EEG_clean.clean_epochs = [clean_epochs_table, clean_epoch_info];

% Calculate the percentage of clean epochs
clean_epochs_percent = (clean_epochs_count / num_epochs) * 100;

% Display the number of clean epochs
disp(['Number of clean epochs: ', num2str(clean_epochs_count)]);
disp(['Percentage of clean epochs: ', num2str(clean_epochs_percent), '%']);
    %% initialize things for spectrum and so on
% stfft sepctrogram                 
    srate = EEG_clean.srate;                 % signal rate
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

    %% loop over channels
    nr_chan = length(EEG_clean.chanlocs);
    if nr_chan > 0
        for chan = 1:nr_chan
            % fulldata = EEG_clean.data(chan,:); % data to go through
            %% funktioniert
            duration_fft = duration_window;    % duration_window or higher length(fulldata)
            
                % check if there are any clean epochs
                nr_epochs = height(EEG_clean.clean_epochs);
                    if nr_epochs > 0
                        for epochs = 1:nr_epochs
                            start_index = EEG_clean.clean_epochs{epochs, 'start_epoch_index'}; %(epochs-1)*epoch_windowlen*(1-epoch_windowover)*EEG_clean.srate + start_index_o; % (epochs-1)*epoch_windowlen*EEG_clean.srate + start_index_o;
                            end_index = EEG_clean.clean_epochs{epochs, 'end_epoch_index'}; % min(start_index + epoch_windowlen*EEG_clean.srate, end_index_o) % min(epochs*epoch_windowlen*EEG_clean.srate + start_index_o, end_index_o);
                            indexs = start_index:end_index;
                            % previously everything from start_index:end_index
                            % indexs = find(ismember(index_cumsum,start_index:end_index) & index_cut); % start_index:end:index
                            %% bands
                            % werte checken "O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_laura_ghb\scripts\ghb_spectrum_jannickmat.m"
                            delta_b = [0.5, 4];
                            theta_b = [4, 8];
                            alpha_b = [8, 12];
                            beta_b = [15, 30];
                            gamma_b = [30, 40];
                            spindles_b = [12, 16];
                           
                            % bands #new
                            delta_ts = bandpass(EEG_clean.data(chan, indexs), delta_b, EEG_clean.srate);
                            theta_ts = bandpass(EEG_clean.data(chan, indexs), theta_b, EEG_clean.srate);
                            alpha_ts = bandpass(EEG_clean.data(chan, indexs), alpha_b, EEG_clean.srate);
                            beta_ts = bandpass(EEG_clean.data(chan, indexs), beta_b, EEG_clean.srate);
                            gamma_ts = bandpass(EEG_clean.data(chan, indexs), gamma_b, EEG_clean.srate);


                            %% bands burst parameters
                            t_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(EEG_clean.data(chan, indexs)), ...
                                fs=EEG_clean.srate, dual_thresh=[1, 1.5], f_range=theta_b);
                            d_t_b = diff([0, double(t_burst.tolist()), 0]);
                            if sum(abs(d_t_b))/2 > 2
                                midpp = (find(d_t_b == 1) + find(d_t_b == -1) - 1) / 2;
                                t_bursts_ssd = std(diff(midpp)/EEG_clean.srate);
                            else
                                t_bursts_ssd = NaN;
                            end
                            t_burst_stat = py.neurodsp.burst.compute_burst_stats(t_burst, fs=EEG_clean.srate);
                            t_burst_stat = struct(t_burst_stat);
                            t_burst_stat.n_bursts = double(t_burst_stat.n_bursts);

                            a_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(EEG_clean.data(chan, indexs)), ...
                                fs=EEG_clean.srate, dual_thresh=[1, 1.5], f_range=alpha_b);
                            d_a_b = diff([0, double(a_burst.tolist()), 0]);
                            if sum(abs(d_a_b))/2 > 2
                                midpp = (find(d_a_b == 1) + find(d_a_b == -1) - 1) / 2;
                                a_bursts_ssd = std(diff(midpp)/EEG_clean.srate);
                            else
                                a_bursts_ssd = NaN;
                            end
                            a_burst_stat = py.neurodsp.burst.compute_burst_stats(a_burst, fs=EEG_clean.srate);
                            a_burst_stat = struct(a_burst_stat);
                            a_burst_stat.n_bursts = double(a_burst_stat.n_bursts);

                            b_burst = py.neurodsp.burst.detect_bursts_dual_threshold(py.numpy.array(EEG_clean.data(chan, indexs)), ...
                                fs=EEG_clean.srate, dual_thresh=[0.9, 1.5], f_range=beta_b);
                            d_b_b = diff([0, double(b_burst.tolist()), 0]);
                            if sum(abs(d_b_b))/2 > 2
                                midpp = (find(d_b_b == 1) + find(d_b_b == -1) - 1) / 2;
                                b_bursts_ssd = std(diff(midpp)/EEG_clean.srate);
                            else
                                b_bursts_ssd = NaN;
                            end
                            b_burst_stat = py.neurodsp.burst.compute_burst_stats(b_burst, fs=EEG_clean.srate);
                            b_burst_stat = struct(b_burst_stat);
                            b_burst_stat.n_bursts = double(b_burst_stat.n_bursts);

                            %% bands coherence within #new
                            t_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=theta_b);
                            a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=alpha_b);
                            t_a_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=[theta_b(1), alpha_b(2)]);
                            b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=beta_b);
                            a_b_lagcoh = py.neurodsp.rhythm.compute_lagged_coherence(py.numpy.array(EEG_clean.data(chan, indexs)), fs=EEG_clean.srate, freqs=[alpha_b(1), beta_b(2)]);


           
                            %% geht nicht hurst long
                            hurs_exp_long = hurst_estimate(zscore(EEG_clean.data(chan, indexs)), 'absval', 0, 1);
                            hurs_alpha = hurst_estimate(zscore(alpha_ts), 'absval', 0, 1);
                            hurs_beta = hurst_estimate(zscore(beta_ts), 'absval', 0, 1);
                            hurs_theta = hurst_estimate(zscore(theta_ts), 'absval', 0, 1);
                            hurs_delta = hurst_estimate(zscore(delta_ts), 'absval', 0, 1);
                             %#new
                            hurs_gamma = hurst_estimate(zscore(gamma_ts), 'absval', 0, 1);

                            %% other fractal dimensions
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
                            
                            %% pwelch for power spectrum
                            % set up pwelch method for power spectrum
                            duration_window = 4*EEG_clean.srate;           % 4s windows
                            which_window = hanning(duration_window); % hanning window
                            overlap_window = duration_window*0.5;    % 50% overlap, 2s

                            % setup epochs
                            epoch_windowlen = 30; % seconds
                            epoch_windowover = 0; % no overlap

                         
                                % do pwelch spectrum
                                [psd, freq] = pwelch(EEG_clean.data(chan,idxs), ...
                                    which_window, overlap_window, EEG_clean.srate/freq_res, EEG_clean.srate);

                                freq_of_interest = freq >= 0.5 & freq <= 30;
                                psdl = 10*log10(psd(freq_of_interest,:))';
                                psdo = psd(freq_of_interest,:)';
                                freql = freq(freq_of_interest,:)';

                                id_swa = find(freql >= delta_b(1) & freql <= delta_b(2));
                                eeg_d_abs = trapz(freql(id_swa), psdo(id_swa));
                                eeg_d_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                                id_swa = find(freql >= theta_b(1) & freql <= theta_b(2));
                                eeg_t_abs = trapz(freql(id_swa), psdo(id_swa));
                                eeg_t_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                                id_swa = find(freql >= alpha_b(1) & freql <= alpha_b(2));
                                eeg_a_abs = trapz(freql(id_swa), psdo(id_swa));
                                eeg_a_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                                id_swa = find(freql >= beta_b(1) & freql <= beta_b(2));
                                eeg_b_abs = trapz(freql(id_swa), psdo(id_swa));
                                eeg_b_rel = trapz(freql(id_swa), psdo(id_swa))/trapz(freql,psdo);

                                %anpassen
                                id_swa = find(freql >= swa_b(1) & freql <= swa_b(2));
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
                                stagecurr = string(EEG_clean.clean_epochs.stage(EEG_clean.clean_epochs.epoch(epochs) == epo));
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
                                    num2cell([eeg_d_abs, eeg_d_rel,...
                                    eeg_t_abs, eeg_t_rel,...
                                    eeg_a_abs, eeg_a_rel,...
                                    eeg_b_abs, eeg_b_rel,...
                                    eeg_s_abs, eeg_s_rel]), ...
                                    num2cell(psdo)], 'VariableNames', colnames_spect);
                            

                        temp_t = temp_t(~cellfun('isempty', temp_t));
                        spectrum_info = [spectrum_info; vertcat(temp_t{:})];

                            
                            %% geht
                            perm_en = PermEn(symb_king(EEG_clean.data(chan, indexs)), 7);
                            perm_en = perm_en(7);
                            perm_en_t = PermEn(symb_king(theta_ts), 7);
                            perm_en_t = perm_en_t(7);
                            perm_en_a = PermEn(symb_king(alpha_ts), 7);
                            perm_en_a = perm_en_a(7);
                            perm_en_b = PermEn(symb_king(beta_ts), 7);
                            perm_en_b = perm_en_b(7);
                            
                            %% geht nicht, googeln?
                            % phase amp coupling #new
                            EEG_cleant = EEG_clean;
                            EEG_cleant.data = EEG_cleant.data(:,indexs);
                            evalc('EEG_cleant = eeg_checkset(EEG_cleant)');
                            EEG_cleant = pop_pac(EEG_cleant, 'channels', theta_b, alpha_b, chan,chan, 'forcecomp', 1); % pop_plotpac(EEG);
                            phase_amp_t_a = mean(EEG_cleant.etc.eegpac(1).glm.pacval, 'all');
                            EEG_cleant.etc.eegpac = [];
                            EEG_cleant.etc.pacplotopt = [];
                            EEG_cleant = pop_pac(EEG_cleant, 'channels', alpha_b, beta_b, chan,chan, 'forcecomp', 1);
                            phase_amp_a_b = mean(EEG_cleant.etc.eegpac(1).glm.pacval, 'all');
                            EEG_cleant.etc.eegpac = [];
                            EEG_cleant.etc.pacplotopt = [];
                            EEG_cleant = pop_pac(EEG_cleant, 'channels', theta_b, beta_b, chan,chan, 'forcecomp', 1);
                            phase_amp_t_b = mean(EEG_cleant.etc.eegpac(1).glm.pacval, 'all');
                            EEG_cleant.etc.eegpac = [];
                            EEG_cleant.etc.pacplotopt = [];
                            clear EEG_cleant;
            
                           
                            %% other time vars
                            kurt = kurtosis(EEG_clean.data(chan, indexs));
                            skew = skewness(EEG_clean.data(chan, indexs));
                            [hjorth_activity, hjorth_mobility, hjorth_complexity] = hjorth(EEG_clean.data(chan, indexs)',0);
            
                            %% other complexity measures
                      
                            [MSx,smse] = MSEn(EEG_clean.data(chan, indexs), MSobject('SampEn')); % multiscale sample enntropy

                            scalz = unique(2*floor(linspace(EEG_clean.srate/30, EEG_clean.srate/0.7, 5)/2)+ 1);
                            mlz = getMultiscaleLZ(EEG_clean.data(chan, indexs), scalz);
                            smlz = sum(mlz);
                            figure();
                            plot(mlz);

                            %% geht nicht
                            % mutual information in bands
                            % paper zu gcmi_cc anschauen
                            % für P channels auch noch machen
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

                                gcmi_cc_helper = @(chan_1, chan_2, remove_channel_which, EEG_clean, indexs) ...
                                    iff(checkzero(remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan_1))) && ...
                                    checkzero(remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan_2))), ...
                                    @() gcmi_cc(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan_1), indexs), ...
                                            EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan_2), indexs)), ...
                                    @() NaN, {}, {});

                                % remove no channels
                                remove_channel_which = zeros(1,8);

                                mi_F_O_within = nanmean([gcmi_cc_helper('O1', 'F3', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('O2', 'F4', remove_channel_which, EEG_clean, indexs)]);
                                mi_F_O_crossed = nanmean([gcmi_cc_helper('O1', 'F4', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('O2', 'F3', remove_channel_which, EEG_clean, indexs)]);
                              
                                mi_F_C_within = nanmean([gcmi_cc_helper('C3', 'F3', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('C4', 'F4', remove_channel_which, EEG_clean, indexs)]);
                                mi_F_C_crossed = nanmean([gcmi_cc_helper('C3', 'F4', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('C4', 'F3', remove_channel_which, EEG_clean, indexs)]);

                                mi_C_O_within = nanmean([gcmi_cc_helper('O1', 'C3', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('O2', 'C4', remove_channel_which, EEG_clean, indexs)]);
                                mi_C_O_crossed = nanmean([gcmi_cc_helper('O1', 'C4', remove_channel_which, EEG_clean, indexs), ...
                                    gcmi_cc_helper('O2', 'C3', remove_channel_which, EEG_clean, indexs)]);


                                mi_F_hemi = gcmi_cc_helper('F3', 'F4', remove_channel_which, EEG_clean, indexs);
                                mi_C_hemi = gcmi_cc_helper('C3', 'C4', remove_channel_which, EEG_clean, indexs);
                                mi_O_hemi = gcmi_cc_helper('O1', 'O2', remove_channel_which, EEG_clean, indexs);

                                
            
                        %% geht wieder
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
            
                             % strongest peak detection in certain band
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
            
                            
                            %% hemispheric spectral measure
                            % wie ähnlich sind zwei spektren (fläche zwischen den spectrograms)
                            if remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'F4'))|remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'F3'))
                                wasserstein_frontal = NaN;
                            else
                                if all(~strcmp({EEG_clean.chanlocs.labels}, 'F4'))|all(~strcmp({EEG_clean.chanlocs.labels}, 'F3'))
                                    wasserstein_frontal = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'F4'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'F3'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_frontal = sum(abs(cf3-cf4))/length(cf3);
                                end
                            end
            
                            if remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'C4'))|remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'C3'))
                                wasserstein_central = NaN;
                            else
                                if all(~strcmp({EEG_clean.chanlocs.labels}, 'C4'))|all(~strcmp({EEG_clean.chanlocs.labels}, 'C3'))
                                    wasserstein_central = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'C4'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'C3'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_central = sum(abs(cf3-cf4))/length(cf3); 
                                end
                            end
            
                            if remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'O2'))|remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, 'O1')) 
                                wasserstein_occipital = NaN;
                            else
                                if all(~strcmp({EEG_clean.chanlocs.labels}, 'O2'))|all(~strcmp({EEG_clean.chanlocs.labels}, 'O1'))
                                    wasserstein_occipital = NaN;
                                else
                                    [psdf4, freqf4] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'O2'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    [psdf3, freqf3] = pwelch(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels},'O1'), indexs), ...
                                                         which_window, overlap_window, duration_fft, srate);
                                    cf3 = cumtrapz(psdf3(freqf3 >= 0.5 & freqf3 <= 30));
                                    cf3 = cf3./cf3(end);
                                    cf4 = cumtrapz(psdf4(freqf4 >= 0.5 & freqf4 <= 30));
                                    cf4 = cf4./cf4(end);
                                    wasserstein_occipital = sum(abs(cf3-cf4))/length(cf3);
                                end
                            end
            
                            
            
            
            
                            % hilbert envelope hurst etc
                            hil = hilbert(EEG_clean.data(chan, indexs));
                            hil_env = abs(hil);
            
                            hurst_hilbert = hurst_estimate(hil_env, 'absval', 0, 1);
                            slope_en_hilbert =  SlopEn(zscore(hil_env));
                            std_hilbert =  std(hil_env);
                            mean_hilbert =  mean(hil_env);
                            rmssd_hilbert = sqrt(mean(diff(hil_env).^2));
                            lzw_hilbert = lzw(hil_env>mean_hilbert);
            
        
        
        %% geht wieder
                            % ECG measures
                            
                            [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(EEG_clean.data_ECG(1,indexs)),srate,0);
                            RR = diff(qrs_i_raw)/srate; % RR in seconds

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
                                py_peaks_result = py.neurokit2.ecg_peaks(py.numpy.array(EEG_clean.data_ECG(1,indexs)), ...
                                    sampling_rate=py.int(srate));
                                py_ecg_rate_result = py.neurokit2.ecg_rate(py_peaks_result, sampling_rate=py.int(srate), desired_length=py.int(length(EEG_clean.data_ECG(1,indexs))));
                                edr_result = double(py.neurokit2.ecg_rsp(py_ecg_rate_result, sampling_rate=srate).tolist());
                                [edr_peak_amp, edr_peak_ind] = findpeaks(edr_result);
                                br_permin = length(edr_peak_ind)/(length(indexs)/srate/60);
                            catch ME
                                br_permin = NaN;
                            end
                               
                            %% save the data in parpool
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
                            %% save data old
                            swa_data = [swa_data; table(string(thename), string(EEG_clean.chanlocs(chan).labels), ...
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






###### again other script, with SW detection and travelling waves
if size(EEG_clean.data,1) > 0 & ~isempty(EEG_clean.scoring)
        nr_epochs = min(ceil((size(EEG_clean.data,2)/EEG_clean.srate)/30),size(current_scoring,1));
        if ~ ceil((size(EEG_clean.data,2)/EEG_clean.srate)/30) == size(current_scoring,1)
            fprintf('------------Warning nr_epochs not matching for iteration %d\n', i);
        end
        nr_channels = size(EEG_clean.data,1);

     
        epoch_windowlen = 30; % seconds
        epoch_windowover = 0; % no overlap
        if length(EEG_clean.data(1, :))/EEG_clean.srate >= 2
            % add stage to each data sample
            current_stages = str2double(string(current_scoring.stage)');
            lastvalue = find(~isnan(table2array(current_scoring(:,'interval_s')))');
            lastvalue = lastvalue(end);
            hypnodata = [repelem(current_stages(1),table2array(current_scoring(1,'start_s'))*EEG_clean.srate), ...
                        repelem(current_stages(1:lastvalue), 30*EEG_clean.srate), ...
                        repelem(current_stages(lastvalue), max(0,(size(EEG_clean.data,2)/EEG_clean.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG_clean.srate))
                        ];
            current_epoch = str2double(string(current_scoring.epoch)');
            epochdata = [repelem(current_epoch(1),table2array(current_scoring(1,'start_s'))*EEG_clean.srate), ...
                        repelem(current_epoch(1:lastvalue), 30*EEG_clean.srate), ...
                        repelem(current_epoch(lastvalue), max(0,(size(EEG_clean.data,2)/EEG_clean.srate- (table2array(current_scoring(lastvalue,'start_s'))+table2array(current_scoring(lastvalue,'interval_s'))))*EEG_clean.srate))
                        ];
                    
                                
            data_for_py = cell(1, size(EEG_clean.data, 1));
            for row = 1:size(EEG_clean.data, 1)
                data_for_py(row) = {EEG_clean.data(row, :)};
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




    %% average pro sleep stage?
    %erster fooof löschen weil zweiter detaillierter?
    % ideen für sleep EEG?
    %yasa spindles?
    % cap pattern (cyclic, acyclic pattern)
    % arousals über ganze nacht, nicht nur clean epochs?
    
