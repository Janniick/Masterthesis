function [bad_channel_idx, value_chan] = eeg_badchannel_detect(EEG, rem_ind_all, cutoff_percent)
arguments
    EEG
    rem_ind_all (1,:) = zeros(1,size(EEG.data,2))
    cutoff_percent (1,1) {mustBeNumeric,mustBeReal} = 0.55
end

% correlations between channels, without vertical lines
corE = abs(corrcoef(EEG.data(:,~rem_ind_all)'));

% a combination of median and sum correlations for each channel 
value_chan = mean(corE);

% percentage of the maximum value
bad_channel_idx = value_chan < quantile(value_chan,0.8)*cutoff_percent;

end