% author : Benjamin Stucky
% year   : 2024
function EEG = eeg_jelena_ecg_removal(EEG, ECG, cut_chunks)
	% EEG: is a eeglab structure with eegchannels, that one wants to remove ecg artifacts
	% ECG: is again an eeglab structure with one ecg channel
	% cut_chunks: is an object with two columns (start, end) in seconds, indicating if there are some strong movement artifacts to exclude from ECG removal, 
	%             rows indicate multiple such artifacts (if occured), this is done, since EEG values at ECG peaks, might be very high and thus the mean to subtract as well, lead to ringing in the EEG.
	% requires pan_tompkin, jelena_ecg_helper, 

 % Jelena's EKG removal technique
    fprintf('---Jelenas ECG removal technique\n');
    % detect R peaks
    if ~isempty(ECG.data)
        [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(ECG.data(1,:)),ECG.srate,0);

        % % delay is not true of pan_thompkin(), do a search oneself
        % delayf = @(x) jelena_ecg_helper(x, qrs_i_raw, ECG); % delayf = @(x) mean(ECG.data(1,qrs_i_raw+x));
        % searchf = -floor(ECG.srate*0.04):floor(ECG.srate*0.04);
        % minf = arrayfun(delayf, searchf);
        % [mf,mif] = max(minf);
        %
        % % actual peak
        % delay_actual = searchf(mif);
        % rpeaks_at = qrs_i_raw + delay_actual;

        % delay is not true of pan_thompkin(), do a search oneself
        delayf = @(x) jelena_ecg_helper(x, qrs_i_raw, ECG); % delayf = @(x) mean(ECG.data(1,qrs_i_raw+x));
        ECGm = ECG;
        ECGm.data = -ECGm.data;
        delayfminus = @(x) jelena_ecg_helper(x, qrs_i_raw, ECGm); % delayf = @(x) mean(ECG.data(1,qrs_i_raw+x));
        searchf = -floor(ECG.srate*0.04):floor(ECG.srate*0.04);

        minf = arrayfun(delayf, searchf);
        [mf,mif] = max(minf);

        minf_m = arrayfun(delayfminus, searchf);
        [mf_m,mif_m] = max(minf_m);

        % actual peak
        delay_actual = searchf(mif);
        rpeaks_at = qrs_i_raw + delay_actual;
        rpeaks_at = rpeaks_at(rpeaks_at>0);
        rpeaks_at = rpeaks_at(rpeaks_at<=size(ECG.data,2));

        delay_actual_m = searchf(mif_m);
        rpeaks_at_m = qrs_i_raw + delay_actual_m;
        rpeaks_at_m = rpeaks_at_m(rpeaks_at_m>0);
        rpeaks_at_m = rpeaks_at_m(rpeaks_at_m<=size(ECGm.data,2));

        [polar_max, polar_idx] = max([mean(ECG.data(1,rpeaks_at)),...
            mean(ECGm.data(1,rpeaks_at_m))]);
        if polar_idx == 2
            disp('ECG polarity switched');
            ECG = ECGm;
            [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(double(ECG.data(1,:)),ECG.srate,0);
            delayf = @(x) jelena_ecg_helper(x, qrs_i_raw, ECG);
            minf = arrayfun(delayf, searchf);
            [mf,mif] = max(minf);
            delay_actual = searchf(mif);
            rpeaks_at = qrs_i_raw + delay_actual;
        end
        clear ECGm;


        % above 400uV
        rpeaks_at = rpeaks_at(abs(ECG.data(1,rpeaks_at))>400);

        % segments around peaks
        segment_len = ECG.srate*1/4; % 1/4 second data around R peak

        %EEG = EEG;
        suborno = zeros(1, size(EEG.data,1));
        for chcha = 1:size(EEG.data,1)
            % Initialize an empty matrix to store EEG segments
            EEG_segments = zeros(length(rpeaks_at), segment_len+1);

            % Extract EEG segments around R peaks, change 1 to channel one likes
            for ip = 1:length(rpeaks_at)
                winlow = (rpeaks_at(ip) - segment_len/2);
                appendlow = [];
                winhigh = (rpeaks_at(ip) + segment_len/2);
                appendhigh = [];
                if(winlow < 1)
                    appendlow = zeros(1,1-winlow);
                    winlow = 1;
                end
                seeg = size(EEG.data,2);
                if(winhigh > seeg)
                    appendhigh = zeros(1,winhigh - seeg);
                    winhigh = seeg;
                end
                EEG_segments(ip, :) = [appendlow, EEG.data(chcha, winlow:winhigh), appendhigh];
            end


            % which peaks to take out
            peaks_to_cut = [];
            if size(cut_chunks{chcha},1)>0
                for chunklen = 1:size(cut_chunks{chcha},1)
                    peaks_to_cut = [peaks_to_cut, find(rpeaks_at/EEG.srate <= cut_chunks{chcha}(chunklen,2) & rpeaks_at/EEG.srate >= cut_chunks{chcha}(chunklen,1))];
                end
            end

            idx_remove = rpeaks_at(peaks_to_cut);

            % moving 30% trimmed mean ECG patterns in the EEG
            avg_EKG_patterns = zeros(size(EEG_segments));
            for ipk = 31:(length(rpeaks_at)-30) % 30 beats to each side
                take_idx = (ipk-30):(ipk+30);
                take_idx = setdiff(take_idx, idx_remove);
                if length(take_idx)>=8 % if more than 8 beats are remaining in an interval go with averaging
                    avg_EKG_patterns(ipk, :) = trimmean(EEG_segments(take_idx, :), 30,1); % mean(EEG_segments((ipk-30):(ipk+30), :), 1);
                end
            end


            % also check if it is cutoff by the band algorithm
            % if size(cut_chunks{chcha},1)>0
            %     for chunklen = 1:size(cut_chunks{chcha},1)
            %         peaks_to_cut = find(rpeaks_at/EEG.srate <= cut_chunks{chcha}(chunklen,2) & rpeaks_at/EEG.srate >= cut_chunks{chcha}(chunklen,1));
            %         peaks_to_cut = max(min(peaks_to_cut)-31,1):min(max(peaks_to_cut)+30,max(rpeaks_at));
            %         avg_EKG_patterns(peaks_to_cut, :) = 0;
            %     end
            % end


            % taper the moving average ECG by Tukey 50%
            avg_EKG_patterns2 = avg_EKG_patterns .* tukeywin(width(avg_EKG_patterns))';

            % figure;
            % ylim([-30 15]);
            % hold on; % Allows multiple lines on the same plot
            % for i = 1:size(avg_EKG_patterns2, 1)
            %     plot((1:size(avg_EKG_patterns2, 2))/EEG.srate, avg_EKG_patterns2(i, :));
            % end
            % hold off;
            %
            % figure;
            % plot((1:(segment_len + 1))/EEG.srate - segment_len/2/EEG.srate ,avg_EKG_patterns2');
            % ylim([-30 15]); % set the limits of the y-axis
            % xlim([-0.125, 0.125]);
            % xlabel('Time (s)')
            % ylabel('EEG pattern (µV)')
            %
            % figure;
            % plot((1:(segment_len + 1))/EEG.srate - segment_len/2/EEG.srate, EEG.data(3,(rpeaks_at(i) - segment_len/2) : (rpeaks_at(i) + segment_len/2)))
            % hold on;
            % plot((1:(segment_len + 1))/EEG.srate - segment_len/2/EEG.srate, EEG.data(3,(rpeaks_at(i) - segment_len/2) : (rpeaks_at(i) + segment_len/2)) - avg_EKG_patterns2(i, :))
            % plot((1:(segment_len + 1))/EEG.srate - segment_len/2/EEG.srate, avg_EKG_patterns2(i, :))
            % xlim([-0.125, 0.125]);
            % xlabel('Time (s)')
            % ylabel('EEG (µV)')
            % legend('Original EEG','Corrected EEG', 'Subtracted EEG')
            % hold off;

            % Subtract the average EKG patterns from the EEG trace at the position of the corresponding R peak
            suborno(chcha) = mad(reshape(avg_EKG_patterns(:,[1:floor((segment_len+1)/3), ...
                2*floor((segment_len+1)/3):3*floor((segment_len+1)/3)]), 1, []))*1.3 < ...
                mad(reshape(avg_EKG_patterns(:,[floor((segment_len+1)/3):2*floor((segment_len+1)/3)]), 1, []));


            suborno(chcha) = 1; % here just run it always to see, othwerwise comment out
            if(suborno(chcha))
                for ips = 1:length(rpeaks_at)
                    winlow = (rpeaks_at(ips) - segment_len/2);
                    takeoutl = 0;
                    winhigh = (rpeaks_at(ips) + segment_len/2);
                    takouth = 0;
                    if(winlow < 1)
                        takeoutl = 1-winlow;
                        winlow = 1;
                    end
                    seeg = size(EEG.data,2);
                    if(winhigh > seeg)
                        takouth = winhigh - seeg;
                        winhigh = seeg;
                    end


                    EEG.data(chcha, winlow:winhigh) = EEG.data(chcha, winlow:winhigh) - avg_EKG_patterns2(ips,(1+takeoutl):(segment_len+1-takouth));
                end
            end
        end
    else
        suborno = zeros(1, size(EEG.data,1));
    end


end