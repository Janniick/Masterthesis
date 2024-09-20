% author : Benjamin Stucky
% year   : 2024
function EEG = eeg_electricalnoise_reducer_notparallel(EEG, omits, lowestfreq,highestfreq, display_plot, plot_output_path, filename, where_notch, type, sequential_min, times_max_mad)
    arguments
        EEG
        omits (1,:) {mustBeNumericOrLogical} = repelem(false, size(EEG.data,2)) % vector size(EEG.data,2) indicating with 1/true where omit the data
        lowestfreq (1,1) {mustBeNumeric} = 5 % below this frequency in Hz no electical noise reduction
        highestfreq (1,1) {mustBeNumeric} = 200 % above this frequency in Hz no electrical noise reduction
        display_plot (1,1) {mustBeNumericOrLogical} = false
        plot_output_path (1,:) {mustBeText} = '' % if '' then no plot output, if string then saved
        filename (1,:) {mustBeText} = 'electricnoisefilter'
        where_notch (1,:) {mustBeNumeric} = [50, 100, 150]
        type (1,:) {mustBeText} = 'median' % type is either median (amplitude to median) zero (setting noise frequencies to 0 + 0i) or bandpass (some narrowband filtering)
        sequential_min (1,1) {mustBeNumeric} = size(EEG.data(:,~omits),2)/EEG.srate/60
        times_max_mad (1,1) {mustBePositive, mustBeNumeric} = 6
    end
    % example
    % eeg_electricalnoise_reducer(EEG, false, '', '', [50, 100, 150], 'median', 15)

    % print what one is doing
    fprintf('---remove electrical noise: Benjis algorithm\n');
    
    % indices to keep
    keep_idx = find(~omits);
    length_keep = length(keep_idx);
    

    eeg_channel_indices = 1:size(EEG.data, 1);

    if display_plot == 1 | ~strcmp(plot_output_path, '')
        if display_plot == 1
            fig = figure('Visible', 'on');
        else
            fig = figure('Visible', 'off');
        end
    end

    for chaneli = eeg_channel_indices
        fprintf('%s', EEG.chanlocs(chaneli).labels);
    
        % notch filter order and zero padding to warm up the filter (10 seconds)
        nf_order = 300;
        pad_length = 10*EEG.srate;
    
        % notch filter
        for notch_run = [1, 0.5, 0.25, 0.1, 0.1] % how much +- Hz
            for notch_freq = where_notch(where_notch<EEG.srate/2)
                notch_filt = designfilt('bandstopfir', ...
                    'FilterOrder', nf_order, ...
                    'CutoffFrequency1', notch_freq-notch_run, ...
                    'CutoffFrequency2', notch_freq+notch_run, ...
                    'SampleRate', EEG.srate);
                filtered_notch = filtfilt(notch_filt, [zeros(1, pad_length), double(EEG.data(chaneli, :))]);
                EEG.data(chaneli,:) = filtered_notch(pad_length + 1:end);
            end
        end
        clear filtered_notch;

        
        nrsequences = max(floor((length_keep/EEG.srate/60)/sequential_min),1);
        %parfor
        for seqi = 1:nrsequences
            % sequence from to index
            fromto = ((seqi-1)*sequential_min*60*EEG.srate+1):seqi*sequential_min*60*EEG.srate;

            % get info on FFT and do FFT
            duration_fft = length(EEG.data(chaneli,keep_idx(fromto)));

            y = fft(EEG.data(chaneli,keep_idx(fromto)));

            magnitude_logy = 10*log10(abs(y(1:floor(duration_fft/2))).^2);
            phase_y = angle(y);

            xw = EEG.srate/duration_fft*(1:floor(duration_fft/2));

            % rolling median
            rolling_window = 2*duration_fft/EEG.srate; % 2Hz steps
            yrollmedian = movmedian(magnitude_logy, rolling_window);

            % cutoff
            % rolling_cutoff = yrollmedian + 4.5*(movmad(magnitude_logy, rolling_window)); % old cutoff
            ymax = movmax(magnitude_logy, (0.1/2)*rolling_window); % 0.1 Hz
            ymax_med = movmedian(ymax, rolling_window);
            rolling_cutoff = ymax_med + times_max_mad*(mad(ymax-ymax_med));

            freq_to_regulate = xw(magnitude_logy > rolling_cutoff);
            freq_to_regulate = freq_to_regulate(freq_to_regulate > lowestfreq);
            freq_to_regulate = freq_to_regulate(freq_to_regulate < highestfreq);

            % no regulation 0.5 around the notches
            for nono = where_notch
                freq_to_regulate = freq_to_regulate(~(freq_to_regulate > nono - 0.5 & freq_to_regulate < nono + 0.5));
            end

            % regulate per group (subsequent < 0.1 Hz)
            freq_to_regulate = sort(freq_to_regulate);

            if length(freq_to_regulate) > 1
                diffreg = diff(freq_to_regulate);
                if diffreg(1)>0.1
                    groupdiff = [1 (1+grp2idx(cumsum(diffreg>0.1)))'];
                else
                    groupdiff = [1 (grp2idx(cumsum(diffreg>0.1)))'];
                end

                %             tgr = tabulate(groupdiff);
                %             groups_too_large = tgr(tgr(:,2)>10, 1);
                lgroups = length(unique(groupdiff));
                if lgroups > 0
                    fprintf('  n_peaks %d; ', lgroups);
                    for ig = 1:lgroups
                        [fmin, fmax] = bounds(freq_to_regulate(groupdiff == ig));


                        % fmin_idx = find(x == fmin);
                        fmin_idx = find((xw - fmin)>=0, 1,'first');
                        fmax_idx = find((xw - fmax)>=0, 1,'first');
                        % +-0.005 Hz drop window
                        fmin_idx_l = find((xw - (fmin - 0.0005))>=0, 1,'first');
                        fmin_idx_r = find((xw - (fmax + 0.0005))>=0, 1,'first');

                        [mm, idxm] = max(magnitude_logy(fmin_idx_l:fmin_idx_r));
                        fprintf(' %.1f', round(xw(fmin_idx_l+(idxm-1)),1));

                        if strcmp(type, 'zero')
                            % filter out, once per group
                            y(fmin_idx_l:fmin_idx_r) = repelem(0, fmin_idx_r-fmin_idx_l + 1);
                            y((duration_fft-fmin_idx_r):(duration_fft-fmin_idx_l)) = repelem(0, fmin_idx_r-fmin_idx_l + 1);
                        elseif strcmp(type, 'median')
                            amplis = sqrt(10.^(magnitude_logy(fmin_idx_l:fmin_idx_r)./10));
                            thanphase = tan(phase_y(fmin_idx_l:fmin_idx_r));
                            amplidiff = sqrt(10.^((magnitude_logy(fmin_idx_l:fmin_idx_r) - yrollmedian(fmin_idx_l:fmin_idx_r))./10));

                            re_part = amplis./(amplidiff.*sqrt(1+tan(phase_y(fmin_idx_l:fmin_idx_r)).^2));
                            im_part = thanphase.*re_part;

                            % filter out, once per group
                            y(fmin_idx_l:fmin_idx_r) = complex(re_part, im_part);
                            y((duration_fft-fmin_idx_r):(duration_fft-fmin_idx_l)) = flip(complex(re_part, im_part));

                        elseif strcmp(type, 'bandpass')
                            if fmax_idx == fmin_idx
                                fmin_idx = fmin_idx - 1;
                                fmax_idx = fmax_idx + 1;
                                fmin_idx_l = fmin_idx_l - 1;
                                fmin_idx_r = fmin_idx_r + 1;
                            end

                            if fmin_idx == fmin_idx_l
                                fmin_idx_l = fmin_idx_l - 1;
                            end

                            if fmax_idx == fmin_idx_r
                                fmin_idx_r = fmin_idx_r + 1;
                            end

                            % filter out, once per group
                            % Design a bandstop filter with valid specification
                            Fp1 = xw(fmin_idx_l); % Passband edge 1
                            Fst1 = xw(fmin_idx); % Stopband edge 1
                            Fst2 = xw(fmax_idx); % Stopband edge 2
                            Fp2 = xw(fmin_idx_r); % Passband edge 2
                            Ap1 = 0.1; % Passband ripple (dB)
                            Ast = max(0.2,max(magnitude_logy(fmin_idx:fmax_idx) - yrollmedian(fmin_idx:fmax_idx))); % Stopband attenuation (dB)
                            Ap2 = 0.1; % Passband ripple (dB)

                            % Design the filter
                            filter_spec = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2', ...
                                Fp1, Fst1, Fst2, Fp2, Ap1, Ast, Ap2, EEG.srate);
                            bandstop_filter = design(filter_spec, 'cheby1');

                            % Apply the bandstop filter to the signal (again zero pad)
                            filtered_band = filter(bandstop_filter, [zeros(1, pad_length), EEG.data(chaneli,:)]);
                            %                 filtered_band = filtfilt(bandstop_filter.sosMatrix, ...
                            %                     bandstop_filter.ScaleValues, ...
                            %                     double([zeros(1, pad_length), filtered_signal]));
                            EEG.data(chaneli,keep_idx(fromto)) = filtered_band(pad_length + 1:end);

                            clear filtered_band;
                        end


                        % progress bar
                        % if ig == 1
                        %     fprintf('  progress: %d%%', round(ig/lgroups*100));
                        % else
                        %     fprintf(repmat('\b', 1, length(sprintf('  progress: %d%%', round((ig-1)/lgroups*100)))));
                        %     fprintf('  progress: %d%%', round(ig/lgroups*100));
                        % end
                    end
                end
            end

            if strcmp(type, 'zero') | strcmp(type, 'median')
                EEG.data(chaneli,keep_idx(fromto)) = ifft(y,'symmetric');
            end

            fprintf('\n');

        end


    
        % plot
        if display_plot == 1 | ~strcmp(plot_output_path, '')
            fy = fft(EEG.data(chaneli,:));
            magnitude_flogy = 10*log10(abs(fy(1:floor(duration_fft/2))).^2);
            clear fy;
    
            subplot(3, ceil(length(eeg_channel_indices) / 3), find(eeg_channel_indices == chaneli));
            plot(xw, magnitude_logy, 'Color', [0, 0.4470, 0.7410]);
            hold on;
            plot(xw(magnitude_logy>rolling_cutoff), magnitude_logy(magnitude_logy>rolling_cutoff), 'o', 'Color', [0, 0.4470, 0.7410], 'MarkerSize', 2);
            plot(xw, magnitude_flogy,'Color', [0.8500, 0.3250, 0.0980]);
            plot(xw,rolling_cutoff);
            title(sprintf('Channel: %s', EEG.chanlocs(chaneli).labels));
            hold off;
        end
    
    end
    
    % if ~strcmp(plot_output_path, '')
    %     print('saving noise reducer image')
    %     sgtitle(filename);
    %     figtot = gcf;
    %     figtot.Position = [0 0 1800 600];
    %     saveas(figtot, fullfile(plot_output_path, sprintf('%s_artefactremoval.png', filename)));
    % end
end