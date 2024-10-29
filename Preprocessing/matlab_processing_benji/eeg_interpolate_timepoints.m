function [EEG, rem_ind_all] = eeg_interpolate_timepoints(EEG, how_many_missing, rem_ind_all)
% for example how_many_missing = 1:2 if only 1 or 2 missing channel at
% timepoints are allowed to be interpolated
    for ipo = how_many_missing % only 1 mising or 2 missing channels
        idx_ipo = sum(rem_ind_all,1) == ipo; % where ipo missing
        if sum(idx_ipo) > 0
            start_ipos = find(diff([0 idx_ipo]) == 1);
            end_ipos = find(diff([idx_ipo 0]) == -1);
            checkipos = start_ipos == end_ipos;
                for ipos = 1:length(start_ipos)
                    if checkipos(ipos) % if only one value
                        if start_ipos(ipos) == 1
                            EEGt = pop_interp(pop_select(EEG, 'point', ...
                                [start_ipos(ipos),start_ipos(ipos)+1]), ...
                                find(rem_ind_all(:,start_ipos(ipos))));
                            EEG.data(:,start_ipos(ipos):end_ipos(ipos)) = EEGt.data(:,1);
                        else 
                            EEGt = pop_interp(pop_select(EEG, 'point', ...
                                [start_ipos(ipos)-1,start_ipos(ipos)]), ...
                                find(rem_ind_all(:,start_ipos(ipos))));
                            EEG.data(:,start_ipos(ipos):end_ipos(ipos)) = EEGt.data(:,2);
                        end
                    else
                        EEGt = pop_interp(pop_select(EEG, 'point',[start_ipos(ipos), end_ipos(ipos)]), ...
                            find(rem_ind_all(:,start_ipos(ipos))));
                        EEG.data(:,start_ipos(ipos):end_ipos(ipos)) = EEGt.data;
                    end
                end
                rem_ind_all(:,idx_ipo) = zeros(size(rem_ind_all,1),sum(idx_ipo));
        end
    end
end