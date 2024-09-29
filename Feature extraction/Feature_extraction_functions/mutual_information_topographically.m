function mutual_information = mutual_information_topographically(EEG_clean, remove_channel_which, indexs, regions, channels)
    % Function to calculate mutual information for all combinations of given regions and channels, including hemispheric comparisons.
    %
    % Inputs:
    %   - EEG_clean: EEG data structure
    %   - remove_channel_which: logical array indicating removed channels
    %   - indexs: indices of interest for the calculation
    %   - regions: a cell array of region labels (e.g., {'F', 'C', 'O', 'P'})
    %   - channels: a cell array of channel labels (e.g., {'F3', 'F4', 'C3', 'C4', 'O1', 'O2', 'P3', 'P4'})
    %
    % Output:
    %   - mutual_information: a struct containing mutual information for within and crossed comparisons, including hemispheric comparisons.
    
    % Initialize mutual information struct
    mutual_information = struct();

    call = @(fun,par) fun(par{:});
    cellget = @(cell, index) cell{index};
    iff = @(test, truefn, falsefn, truePar, falsePar) call( ...
        cellget({falsefn, truefn}, test + 1), ... % functions
        cellget({falsePar, truePar}, test + 1) ... % params
    );

    checkzero = @(xy) iff(isempty(xy), @()false, @()xy == 0, {},{});

    % Define helper function to calculate mutual information
    gcmi_cc_helper = @(chan_1, chan_2, remove_channel_which, EEG_clean, indexs) ...
        iff(checkzero(remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan_1))) && ...
        checkzero(remove_channel_which(strcmp({EEG_clean.chanlocs.labels}, chan_2))), ...
        @() gcmi_cc(EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan_1), indexs), ...
                EEG_clean.data(strcmp({EEG_clean.chanlocs.labels}, chan_2), indexs)), ...
        @() NaN, {}, {});

    % List of region pairs to compare (within and crossed)
    for i = 1:length(regions)
        for j = i+1:length(regions)
            % Get corresponding channels
            region_1_channels = channels{i};  % e.g., {'F3', 'F4'}
            region_2_channels = channels{j};  % e.g., {'C3', 'C4'}
            
            % Mutual information within (same hemisphere)
            mutual_information.([regions{i}, '_', regions{j}, '_within']) = nanmean([ ...
                gcmi_cc_helper(region_1_channels{1}, region_2_channels{1}, remove_channel_which, EEG_clean, indexs), ...
                gcmi_cc_helper(region_1_channels{2}, region_2_channels{2}, remove_channel_which, EEG_clean, indexs)]);
            
            % Mutual information crossed (opposite hemispheres)
            mutual_information.([regions{i}, '_', regions{j}, '_crossed']) = nanmean([ ...
                gcmi_cc_helper(region_1_channels{1}, region_2_channels{2}, remove_channel_which, EEG_clean, indexs), ...
                gcmi_cc_helper(region_1_channels{2}, region_2_channels{1}, remove_channel_which, EEG_clean, indexs)]);
        end
    end

    % Hemispheric comparisons
    mutual_information.F_hemi = gcmi_cc_helper('F3', 'F4', remove_channel_which, EEG_clean, indexs);
    mutual_information.C_hemi = gcmi_cc_helper('C3', 'C4', remove_channel_which, EEG_clean, indexs);
    mutual_information.O_hemi = gcmi_cc_helper('O1', 'O2', remove_channel_which, EEG_clean, indexs);
    mutual_information.P_hemi = gcmi_cc_helper('P3', 'P4', remove_channel_which, EEG_clean, indexs);

end
