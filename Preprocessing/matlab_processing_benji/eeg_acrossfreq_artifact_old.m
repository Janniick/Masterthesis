% author : Benjamin Stucky
% year   : 2024
function [remove_channel_which, cut_chunks] = eeg_acrossfreq_artifact(EEG, prob_select, percentage_mad_scale, factor_mad_scale_to_moving, moving_window_s)
    arguments
        EEG
        prob_select (1,1) {mustBeNumeric,mustBeReal, mustBeLessThanOrEqual(prob_select,1), mustBeGreaterThanOrEqual(prob_select,0)} = 0.2
		percentage_mad_scale (1,1) {mustBeNumeric,mustBeReal, mustBeLessThanOrEqual(percentage_mad_scale,1), mustBeGreaterThanOrEqual(percentage_mad_scale,0)} = 0.45
        factor_mad_scale_to_moving (1,1) {mustBeNumeric,mustBeReal} = 1
        moving_window_s (1,1) {mustBeInteger,mustBeReal, mustBeGreaterThanOrEqual(moving_window_s,1)} = 240
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove the movement blocks that peak out
    fprintf('---detect blocks to remove from calculation\n');

    % set parameters
    srate = EEG.srate;
    window_length = 1*srate; % 2 sec previously
    noverlap = floor(window_length*0.5);
    padding_sec = (window_length/srate)/2;
    nfft = []; % Use default NFFT
    

    % horizontal line parameters
    movfft_seconds = moving_window_s;
    mad_scale = 2:0.1:4;     % mad_scale = 2.9;  % 2-2.5 for 2s and 0.5 overlap, 3.5;
    mad_scale_mov = factor_mad_scale_to_moving*mad_scale;

    cut_chunks = cell(1, size(EEG.data,1));
    outlier_channel = [];
    fft_per_channel = cell(size(EEG.data,1), 1);
    for chana = 1:size(EEG.data,1)
        % calculate fft, moving functions, to create cutoff
        [sblocks, fblocks, tblocks] = spectrogram(EEG.data(chana,:), hamming(window_length), noverlap, nfft, srate);

        % store it for later bad channel detection type 2
        fft_per_channel{chana, 1} = sblocks(fblocks >= 0.5 & fblocks <= 40, :);

        % for each frequency version
        freq_window_idx = find(fblocks >= 0.5 & fblocks <= 40);
        freqcut = zeros(length(freq_window_idx), size(sblocks,2));
        for fse = 1:length(freq_window_idx)
            % overall version    
            blocks_mag = 10*log10(abs(sblocks(freq_window_idx(fse),:)).^2);
            cutfft_high = median(blocks_mag) + mad_scale*mad(blocks_mag);
            cutfft_low = median(blocks_mag) - mad_scale*mad(blocks_mag)*2;

            movmedfft = movmedian(blocks_mag, movfft_seconds/((window_length/EEG.srate)/(noverlap/EEG.srate)));
            movmadfft = movmad(blocks_mag, movfft_seconds/((window_length/EEG.srate)/(noverlap/EEG.srate)));
            movcutfft_high = movmedfft + mad_scale_mov' .*movmadfft;
            movcutfft_low = movmedfft - mad_scale_mov' .*movmadfft*2;
            finalcutfft_high = min([movcutfft_high,cutfft_high'], [], 2);
            finalcutfft_low = max([movcutfft_low,cutfft_low'], [], 2);

            finalcutt = (finalcutfft_high < blocks_mag)|(finalcutfft_low > blocks_mag);
            
            freqcut(freq_window_idx(fse),:) = sum(finalcutt,1)/size(finalcutt,1)> percentage_mad_scale;
        end
        prob_cut = sum(freqcut,1)/size(freqcut, 1);
        overall_cut = prob_cut > prob_select;

        %         figure();
        %         subplot(2,1,1)
        %         spectrogram(EEG.data(chana,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
        %         ylim([0.5,40])
        %         caxis([-40 40]);
        %         colormap('jet');
        %         colorbar;
        %         subplot(2,1,2)
        %         plot(tblocks/60, prob_cut)

        %         % overall version
        %         blocks_mag = mean(10*log10(abs(sblocks).^2),1);
        %         cutfft = median(blocks_mag) + mad_scale*mad(blocks_mag);
        %         movmedfft = movmedian(blocks_mag, movfft_seconds/((window_length/EEG.srate)/(noverlap/EEG.srate)));
        %         movmadfft = movmad(blocks_mag, movfft_seconds/((window_length/EEG.srate)/(noverlap/EEG.srate)));
        %         movcutfft = movmedfft + mad_scale_mov*movmadfft;
        %         finalcutfft = min(movcutfft, cutfft);
        %         overall_cut = finalcutfft < blocks_mag;
        %
        %         figure();
        %         subplot(2,1,1)
        %         spectrogram(EEG.data(chana,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
        %         ylim([0.5,40])
        %         caxis([-40 40]);
        %         colormap('jet');
        %         colorbar;
        %         subplot(2,1,2)
        %         plot(tblocks/60, mean(10*log10(abs(sblocks).^2),1))
        %         hold on;
        %         plot(tblocks/60, movmedfft);
        %         plot(tblocks/60, finalcutfft);
        %         plot(tblocks([1,end])/60, [cutfft, cutfft]);
        %         hold off;

        % get the cut off chunks start and end in seconds
        wherefftcut = tblocks(overall_cut);
        if length(wherefftcut)>0
            jumpsfft = diff(wherefftcut) > ((window_length/EEG.srate)/(noverlap/EEG.srate));
            wherejump = find(jumpsfft);
            jumpend = wherefftcut([wherejump, end]);
            jumpstart = wherefftcut([1, wherejump+1]);
            cut_chunks{chana} = [max(jumpstart-padding_sec,0); min(jumpend+padding_sec, max(tblocks)+padding_sec)]';

            % per chunk set to NA
            %             xvec = 1:length(EEG.data(chana,:))/EEG.srate;
            %             for chunki = 1:size(cut_chunks,1)
            %                 EEG.data(chana,xvec>=cut_chunks(chunki,1) & xvec <= cut_chunks(chunki,2)) = NaN;
            %             end
        else
            cut_chunks{chana} = [];
        end

        % the stuff that remains, check if it is much larger than the other
        % channels
        tmp_chunk = cut_chunks{chana};
        indices_toremove = [];
        for ikk = 1:size(tmp_chunk, 1)
            indices_toremove = [indices_toremove; find((tblocks-padding_sec) >= tmp_chunk(ikk, 1) & (tblocks+padding_sec) <= tmp_chunk(ikk, 2))']; % padding sec?
        end

        outlier_channel(:,chana) = mean(10*log10(abs(sblocks(fblocks >= 0.5 & fblocks <= 40, setdiff(1:length(tblocks), indices_toremove))).^2),2);

    end


    % !!!!!!!!!!!!!! this changed
    % remove channel!
    % remove_channel_which = sum(outlier_channel > median(outlier_channel,2) + 2.5*mad(outlier_channel,0,2),1)/size(outlier_channel,1) > 0.7 | ...
    % sum(outlier_channel < median(outlier_channel,2) - 4*mad(outlier_channel,0,2),1)/size(outlier_channel,1) > 0.7; % 100* mad??

    sum_outlier = zeros(1,size(outlier_channel,2));
    madiseq = 1:0.1:4;
    for madi = madiseq
        sum_outlier = sum_outlier + sum(outlier_channel > median(outlier_channel,2) + madi*mad(outlier_channel,0,2),1)/size(outlier_channel,1) + ...
            sum(outlier_channel < median(outlier_channel,2) - 1.5*madi*mad(outlier_channel,0,2),1)/size(outlier_channel,1);
    end
    remove_channel_which = sum_outlier/length(madiseq) > 0.7;


    %%%%%%%%%%%%%%%%%%%%%%%%%
    % remove channel, type 2!
    remove_ind_all = [];
    for i_f = 1:size(cut_chunks,2)
        if ~isempty(cut_chunks{i_f})
            for i_i = 1:size(cut_chunks{i_f},1)
                remove_ind_all = [remove_ind_all, (cut_chunks{i_f}(i_i,1)*EEG.srate):(cut_chunks{i_f}(i_i,2)*EEG.srate)];
            end
        end
    end
    remove_ind_all = unique(remove_ind_all);

    matching_windows = [];
    total_windows = floor((size(EEG.data,2) - window_length) / (window_length - noverlap)) + 1;
    start_indices = ((1:total_windows) - 1) * (window_length - noverlap) + 1;
    end_indices = start_indices + window_length - 1;
    for rm_l = remove_ind_all
        matching_windows = [matching_windows, find(rm_l >= start_indices & rm_l <= end_indices)];
    end

    matching_windows = unique(matching_windows);
    only_no_movement_windows = setdiff(1:total_windows, matching_windows);


    fbcut = fblocks(fblocks >= 0.5 & fblocks <= 40);
    which_range =  fbcut >= 0.5 & fbcut <= 25;
    channel_diffs = zeros(size(EEG.data,1),size(EEG.data,1),sum(which_range));
    for i_c = 1:size(fft_per_channel,1)
        for j_c = 1:size(fft_per_channel,1)
            tmp_c = (10*log10(abs(fft_per_channel{i_c,1}(which_range,only_no_movement_windows)).^2) - ...
                10*log10(abs(fft_per_channel{j_c,1}(which_range,only_no_movement_windows)).^2));
            channel_diffs(i_c, j_c,:) = mean(tmp_c,2); % abs(mean(tmp_c(:)))/std(abs(tmp_c(:)))>0.1; %  what to do here!!!
        end
    end

    sum_diffs = zeros(size(channel_diffs));
    ittseq = 0:0.05:10;
    for itt = ittseq
        sum_diffs = sum_diffs + (channel_diffs > itt)*(2-exp(-(0.5*itt).^2));
    end

    helperfun_multiply = @(x) sum(sum_diffs(:,:,x)*sum_diffs(:,:,x)', 1)/(length(ittseq)^2); % this is the core function of interest per frequency x
    sum_multi = cell2mat(arrayfun(helperfun_multiply, 1:size(sum_diffs, 3), 'UniformOutput', false)'); % applies the core function over frequency
    % check all sum_multi if > than a cutoff range @gt means greater than function
    cutoff_range = 1:0.01:20;
    per_freq_final = sum(bsxfun(@gt, sum_multi, reshape(cutoff_range, [1 1 numel(cutoff_range)])), 3) / numel(cutoff_range); % per freq final is percentage of outlier per frequency

    remove_channel_type2 = sum(per_freq_final,1)/size(per_freq_final,1) > 0.35;

    % per sub, old version
    % fbcut = fblocks(fblocks >= 0.5 & fblocks <= 40);
    % channel_diffs = zeros(size(EEG.data,1),size(EEG.data,1));
    % for i_c = 1:size(fft_per_channel,1)
    %     for j_c = 1:size(fft_per_channel,1)
    %         tmp_c = (10*log10(abs(fft_per_channel{i_c,1}(fbcut >= 10 & fbcut <= 40,only_no_movement_windows)).^2) - ...
    %             10*log10(abs(fft_per_channel{j_c,1}(fbcut >= 10 & fbcut <= 40,only_no_movement_windows)).^2));
    %         channel_diffs(i_c, j_c) = mean(tmp_c(:)); % abs(mean(tmp_c(:)))/std(abs(tmp_c(:)))>0.1; %  what to do here!!!
    %     end
    % end
    %
    % sum_diffs = zeros(size(channel_diffs));
    % ittseq = 0:0.05:10;
    % for itt = ittseq
    %     sum_diffs = sum_diffs + (channel_diffs > itt)*(2-exp(-(0.5*itt).^2));
    % end
    %
    % % option 2
    % remove_channel_type3 = sum(sum_diffs*sum_diffs',1)/(length(ittseq)^2) > 6.4; % > 2.2   % or old version % sum(sum_diffs,2)/(length(ittseq)) % > 1
    %
    %
    % remove_channel_which = remove_channel_which | (remove_channel_type2 & remove_channel_type3);

    remove_channel_which = remove_channel_which | remove_channel_type2;

    % hard code 2 falsely detected artifact channels
    if i == 56
        remove_channel_which(4) = 0;
    end
    if i == 88
        remove_channel_which(5) = 0;
    end


    %noisyout = findNoisyChannels(EEG);

    %%%%%%%%%%%%%%%%%%%%%%%%%


    if(sum(remove_channel_which)>0)
        for cit = find(remove_channel_which)
            cut_chunks{cit} = [0, max(tblocks) + padding_sec];
        end
    end

    %%%%% here I could also do merge chunks, if say <4 seconds apart?
    % merge dense chunks together
    cutoff_block = 0.8;
    if size(EEG.data,1) > 0
        for chuju = 1:size(EEG.data,1)
            chunks_temp = cut_chunks{chuju};

            if size(chunks_temp,1) > 0
                chunk_len = chunks_temp(:,2) - chunks_temp(:,1);
                chunk_dist = chunks_temp(2:end,1) - chunks_temp(1:end-1,2);
                full_len = floor((size(EEG.data,2)/srate));
                [~, chunk_sorted] = sort(chunk_dist, 'descend');
                [~, chunk_rank] = sort(chunk_sorted);
                if sum(chunk_len)/full_len > cutoff_block % if total is more black than cutoff percentage remove
                    cut_chunks{chuju} = [0, max(tblocks) + padding_sec];
                else
                    if length(chunk_len) > 1
                        for splitchunk = 1:length(chunk_len)-1
                            splits = cumsum([0;chunk_rank <= splitchunk])+ 1;
                            fsplits = find(accumarray(splits, chunk_len)./(accumarray(splits, chunks_temp(:,2), [], @max) - ...
                                accumarray(splits, chunks_temp(:,1), [], @min)) > cutoff_block);

                            if length(fsplits) > 0
                                for fsi = 1:length(fsplits)
                                    chuchu_find = chunks_temp(splits == fsplits(fsi),:);
                                    cut_chunks{chuju} = unique([cut_chunks{chuju}; min(chuchu_find, [],'all'), max(chuchu_find, [],'all')], 'rows');
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end