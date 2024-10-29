function [EEG, fbch, bads, cut_chunks, rem_ind_all] = eeg_badchannel_interpolate(EEG, bad_chans, cut_chunks, rem_ind_all, channel_names, channel_neighbour_matrix)
% fbch = indices which were interpolated
% bads = indices where no interpolation due to channel_neighbour_matrix
% could be made (needs to have enough neighbours)
if sum(bad_chans)>0
    % make a check here, which of the bad chans can be interpolated
    searchi = ismember(channel_names, {EEG.chanlocs(bad_chans).labels});
    fbch = find(bad_chans);

    available_chans = sum(channel_neighbour_matrix(searchi,~searchi),2)';
    % those cannot be saved, remove all
    bads = fbch(available_chans<=1);
    if length(bads)>0
        for ciu = 1:length(bads)
            cut_chunks{ciu} = [1, length(EEG.data)];
            rem_ind_all(ciu,:) = 1;
        end
    end

    % interpolate those
    fbch = fbch(available_chans>1);
    EEG = pop_interp(EEG, fbch);

    toprint = {EEG.chanlocs(fbch).labels};
    fprintf( 'Interpolated: ');fprintf( '%s ', toprint{:} );fprintf('\n');
else
    fbch = [];
    bads = [];
end

end