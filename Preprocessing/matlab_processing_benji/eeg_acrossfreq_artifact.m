% author : Benjamin Stucky
% year   : 2024
function [remove_channel_which, cut_chunks, rem_ind_all] = eeg_acrossfreq_artifact(EEG, prob_select, percentage_mad_scale, factor_mad_scale_to_moving, moving_window_s, remove_channels, prob_remove_ch)
arguments
    EEG
    prob_select (1,1) {mustBeNumeric,mustBeReal, mustBeLessThanOrEqual(prob_select,1), mustBeGreaterThanOrEqual(prob_select,0)} = 0.2
    percentage_mad_scale (1,1) {mustBeNumeric,mustBeReal, mustBeLessThanOrEqual(percentage_mad_scale,1), mustBeGreaterThanOrEqual(percentage_mad_scale,0)} = 0.45
    factor_mad_scale_to_moving (1,1) {mustBeNumeric,mustBeReal} = 1
    moving_window_s (1,1) {mustBeInteger,mustBeReal, mustBeGreaterThanOrEqual(moving_window_s,1)} = 240
    remove_channels (1,1) {mustBeNumericOrLogical} = false
    prob_remove_ch (1,1) {mustBeNumeric,mustBeReal, mustBeLessThanOrEqual(prob_remove_ch,1), mustBeGreaterThanOrEqual(prob_remove_ch,0)} = 0.3
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the movement blocks that peak out
fprintf('---detect blocks to remove from calculation\n');

% set parameters
srate = EEG.srate;
window_length = 2*srate; % 2 sec previously
noverlap = 0;
padding_sec = 2;
nfft = []; % Use default NFFT


% horizontal line parameters
movfft_seconds = moving_window_s;
mad_scale = 2:0.1:3.5;     % mad_scale = 2.9;  % 2-2.5 for 2s and 0.5 overlap, 3.5;
mad_scale_mov = factor_mad_scale_to_moving*mad_scale;

cut_chunks = cell(1, size(EEG.data,1));
outlierinfo = zeros(size(EEG.data,1), floor(size(EEG.data,2)/window_length));
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

        movmedfft = movmedian(blocks_mag, movfft_seconds);

        % do it for a range of cutoffs, not just a single one

        cutfft_high = repelem(movmedfft,length(mad_scale),1) + (mad_scale*mad(blocks_mag-movmedfft))';
        cutfft_low = repelem(movmedfft,length(mad_scale),1) - (mad_scale*mad(blocks_mag-movmedfft))';


        finalcutt = (cutfft_high < blocks_mag)|(cutfft_low > blocks_mag);

        freqcut(fse,:) = sum(finalcutt,1)/size(finalcutt,1) > percentage_mad_scale;
    end
    prob_cut = sum(freqcut,1)/size(freqcut, 1);
    overall_cut = prob_cut > prob_select;

    outlierinfo(chana, :) = overall_cut;

    % get the cut off chunks start and end in seconds
    wherefftcut = tblocks(overall_cut);
    if length(wherefftcut)>0
        jumpsfft = diff(wherefftcut) > (window_length/EEG.srate);
        wherejump = find(jumpsfft);
        jumpend = wherefftcut([wherejump, end]);
        jumpstart = wherefftcut([1, wherejump+1]);
        cut_chunks{chana} = [max(jumpstart-padding_sec,0); min(jumpend+padding_sec, max(tblocks)+padding_sec)]';
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

    if remove_channels
        outlier_channel(:,chana) = mean(10*log10(abs(sblocks(fblocks >= 0.5 & fblocks <= 40, setdiff(1:length(tblocks), indices_toremove))).^2),2);
    end

end


% !!!!!!!!!!!!!! this changed
% remove channel!
% remove_channel_which = sum(outlier_channel > median(outlier_channel,2) + 2.5*mad(outlier_channel,0,2),1)/size(outlier_channel,1) > 0.7 | ...
% sum(outlier_channel < median(outlier_channel,2) - 4*mad(outlier_channel,0,2),1)/size(outlier_channel,1) > 0.7; % 100* mad??



%%%%%%%%%%%%%%%%%%%%%%%%%
% remove channel, type 2!
rem_wins = outlierinfo;%sum(outlierinfo,1) > 0; 

% indices for whole recording
rem_ind_all =  repelem(rem_wins, 1, window_length); %repelem(rem_wins,window_length);
if size(EEG.data,2)> size(rem_ind_all,2)
    rem_ind_all = [rem_ind_all, repelem(0,size(rem_ind_all,1), size(EEG.data,2)-length(rem_ind_all))]; % remove_ind_all = find(rem_ind_all)
end

% add padding
drem = diff([repelem(0, size(rem_ind_all,1), 1), rem_ind_all, repelem(0, size(rem_ind_all,1), 1)], 1, 2);
for dremi = 1:size(drem,1)
    starti = find(drem(dremi,:) == 1) + 1-floor(padding_sec*EEG.srate);
    endi = find(drem(dremi,:) == -1) -2 + floor(padding_sec*EEG.srate);

    ind_ones = arrayfun(@(s, e) s:e, starti, endi, 'UniformOutput', false);
    ind_ones = [ind_ones{:}];
    ind_ones = ind_ones(ind_ones > 0 & ind_ones <= length(rem_ind_all));
    rem_ind_all(dremi, ind_ones) = 1;
end




% add padding of one window (if wished)
rem_wins = sum(outlierinfo,1) > 0;
drem = diff([0, rem_wins, 0]);
starti = find(drem == 1);
endi = find(drem == -1) - 1;
starti = setdiff(starti,1);
endi = setdiff(endi,length(rem_wins));
rem_wins(starti - 1) = 1;
rem_wins(endi + 1) = 1;

keep_wins = ~rem_wins; % only_no_movement_windows = find(keep_wins)


if remove_channels
    sum_outlier = zeros(1,size(outlier_channel,2));
    madiseq = 1:0.1:4;
    for madi = madiseq
        sum_outlier = sum_outlier + sum(outlier_channel > median(outlier_channel,2) + madi*mad(outlier_channel,0,2),1)/size(outlier_channel,1) + ...
            sum(outlier_channel < median(outlier_channel,2) - 1.5*madi*mad(outlier_channel,0,2),1)/size(outlier_channel,1);
    end
    remove_channel_which = sum_outlier/length(madiseq) > 0.7;


    fbcut = fblocks(fblocks >= 0.5 & fblocks <= 40);
    which_range =  fbcut >= 0.5 & fbcut <= 25;
    u_fft_per_channel = permute(cat(3, fft_per_channel{:}), [3 1 2]);
    u_fft_per_channel = u_fft_per_channel(:,which_range,keep_wins);
    u_fft_per_channel = 10*log10(abs(u_fft_per_channel).^2);
    u_fft_per_channel_expanded = reshape(u_fft_per_channel, [size(u_fft_per_channel,1),1,size(u_fft_per_channel,2:3)]);
    channel_diffs = bsxfun(@minus, u_fft_per_channel_expanded, permute(u_fft_per_channel_expanded, [2, 1, 3, 4]));
    channel_diffs = mean(channel_diffs, 4);



    sum_diffs = zeros(size(channel_diffs));
    ittseq = 0:0.05:10;
    for itt = ittseq
        sum_diffs = sum_diffs + (channel_diffs > itt | channel_diffs < -itt)*(2-exp(-(0.5*itt).^2))/(2*length(ittseq));
    end

    % here check what to do!!!!
    helperfun_multiply = @(x, sumd) sum(sumd(:,:,x)*sumd(:,:,x)'/norm(sumd(:,:,x))/size(sumd, 3), 1); % this is the core function of interest per frequency x
    sum_multi = cell2mat(arrayfun(@(x)helperfun_multiply(x, sum_diffs), 1:size(sum_diffs, 3), 'UniformOutput', false)'); % applies the core function over frequency
    % check all sum_multi if > than a cutoff range @gt means greater than function


    % insert functionality

    % anonymous if function hack
    call = @(fun,par) fun(par{:});
    cellget = @(cell, index) cell{index};
    iff = @(test, truefn, falsefn, truePar, falsePar) call( ...
        cellget({falsefn, truefn}, test + 1), ... % functions
        cellget({falsePar, truePar}, test + 1) ... % params
        );
    inserth = @(a, x, n, idxxx) iff(n==1, @()[a, x(idxxx)], @()iff(n==length(idxxx), @()[x(idxxx), a], @()[x(idxxx(1:(n-1))), a, x(idxxx(n:end))], {},{}), {},{});


    % permutations keep diagonal fixed
    storep = zeros([size(sum_multi),1000]);
    for iit = 1:1000
        perm_diffs = sum_diffs;
        for pi = 1:size(sum_diffs,3)
            idxx = randperm(size(sum_diffs,1)-1);
            withoutii = arrayfun(@(x)setxor(1:size(sum_diffs,1),x), 1:size(sum_diffs,2), 'UniformOutput', false);
            perm_diffs(:,:,pi) = cell2mat(arrayfun(@(x)sum_diffs(inserth(x,withoutii{x},x, idxx),x), 1:length(withoutii), 'UniformOutput', false));
        end
        multip = cell2mat(arrayfun(@(x)helperfun_multiply(x, perm_diffs), 1:size(perm_diffs, 3), 'UniformOutput', false)');
        storep(:,:,iit) = multip;
    end
    per_freq_final = sum_multi > quantile(storep,0.999,3); % mean(quantile(max(storep,[],3), 0.99))

    % cutoff_range = 1:0.01:20;
    % per_freq_final = sum(bsxfun(@gt, sum_multi, reshape(cutoff_range, [1 1 numel(cutoff_range)])), 3) / numel(cutoff_range); % per freq final is percentage of outlier per frequency

    remove_channel_type2 = mean(per_freq_final,1) >= prob_remove_ch;

    remove_channel_which = remove_channel_which | remove_channel_type2;


    %%%%%%%%%%%%%%%%%%%%%%%%%


    if(sum(remove_channel_which)>0)
        for cit = find(remove_channel_which)
            cut_chunks{cit} = [0, max(tblocks) + padding_sec];
        end
    end


else
    remove_channel_which = zeros(1,size(EEG.data,1));
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

                        if ~isempty(fsplits)
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