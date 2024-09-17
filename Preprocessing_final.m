%% All these functions need to be loaded first
% other functions (maybe here not needed, but for feature extraction at least they are needed
addpath(genpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_mesmart\scripts'));
% preprocessing Funktionen
addpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\matlab_processing_benji\');


%%% Add the path to your EEGLAB folder
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui; % Initialize EEGLAB (eeglab nogui; if you want no GUI)

%%% Add the path to your .edf files (if will find all .edf files in all
%%% subfolders)
folderpath = 'D:\Masterarbeit Jannick\Data\GHB_TRA\GHB_TRA_EEG';
edf_files = dir(fullfile(folderpath, '*\*.edf'));
web_files = dir(fullfile(folderpath, '*\*.web'));

%%% define study name (avoid "/" and "_")
study = 'GHBTRA';
%%% define which channels to look at
selected_channels = {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'};
%%% define ASR window length in minutes (default = 8)
windowDuration = 8; 
%%% define ASR SD-Cutoff (default = 30)
SDCutoff = 30;


%% loop through all files
for i = 1:numel(edf_files) %randperm(numel(edf_files)) % 1:numel(edf_files)
    %% import 
    fprintf('********************************************%s****************************************\n', repmat('*', 1, length(edf_files(i).name)));
    fprintf('---------------------------------------File:%s----------------------------------------\n', edf_files(i).name);
    fprintf('                                  Iteration: %s\n', string(i));
    fprintf('                             Percentage run: %s\n', strcat(string(round(100*i/numel(edf_files),0)), '%'));
    fprintf('********************************************%s****************************************\n', repmat('*', 1, length(edf_files(i).name)));
    
    % which filepath
    edf_name = edf_files(i).name;
    edf_path = fullfile(edf_files(i).folder, edf_name);
    % web_name = web_files(i).name;
    % web_path = fullfile(web_files(i).folder, web_name);
    thename = split(edf_name, '.');
    thename = thename{1};
    % Use fileparts to extract the filename without the extension
    [~, name, ~] = fileparts(edf_name);
    % Where to save the generated files
    [pathtosave, ~, ~] = fileparts(edf_path);
    % participant number
    % Get the participant number
    [parentPath, ~, ~] = fileparts(edf_path);
    [~, participant, ~] = fileparts(parentPath);

    % Use eeglab biosig to load edf eeg file
    fprintf('---load data\n');
    EEG_raw = pop_biosig(edf_path);
    EEG_raw.name=edf_name;

    % select channels
    channelLabels = {EEG_raw.chanlocs.labels};
    % Define the pairs of mastoids to check
    pair1 = {'A1', 'A2'};
    pair2 = {'M1', 'M2'};
    % Check if either pair1 or pair2 is present in channelLabels
    isPair1Present = all(ismember(pair1, channelLabels));
    isPair2Present = all(ismember(pair2, channelLabels));
    % Select the channels based on which pair is present
    if isPair1Present
        % If A1 and A2 are present, select them
        EEG_mastoids = pop_select(EEG_raw, 'channel', pair1);
    elseif isPair2Present
        % If M1 and M2 are present, select them
        EEG_mastoids = pop_select(EEG_raw, 'channel', pair2);
    else
        error('Neither A1/A2 nor M1/M2 channels are present in the data.');
    end
    ECG = pop_select(EEG_raw, 'channel', {'ECG'});
    EEG_raw = pop_select(EEG_raw, 'channel', selected_channels);
    newchannelLabels = {EEG_raw.chanlocs.labels};

    % add channel location information
    % on server: readlocs('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_mesmart\eeglab_current\Standard-10-5-Cap385.sfp')
    newlocs = readlocs('D:\Masterarbeit Jannick\scripts\2_Preprocessing\Standard-10-5-Cap385.sfp'); 
    newlocl = {newlocs.labels};
    chanlocsnew = struct();
    fieldsn = fieldnames(newlocs);
    for k = 1:numel(fieldsn)
        chanlocsnew.(fieldsn{k}) = [];
    end
    for f = 1:length({EEG_raw.chanlocs.labels})
        whichrow = find(ismember(newlocl, EEG_raw.chanlocs(f).labels));
        chanlocsnew(f) = newlocs(whichrow);
    end
    EEG_raw.chanlocs = chanlocsnew;


    % downsample
    fprintf('---downsample\n');
    EEG_raw = pop_resample(EEG_raw, 256);
    ECG = pop_resample(ECG, 256);
    EEG_mastoids = pop_resample(EEG_mastoids, 256);
    
    % Save ECG file
    filePath = fullfile(pathtosave, [study, '_', participant, '_ECG.mat']);
    % Save the ECG variable in the specified path
    save(filePath, 'ECG');
    % Optionally display the path where the file was saved
    disp(['Saved ECG in ', filePath]);

    % define parameters for plotting
    screenSize = get(0, 'ScreenSize');
    srate = EEG_raw.srate;
    window_length = 4*srate;
    noverlap = window_length/2;
    nfft = []; % Use default NFFT
    fpass = [0.5 40];
    % how to divide seconds to get values in the plotting frame
    if size(EEG_raw.data(1,:),2) < srate*60 % less than a minute (in seconds)
        divideby = 1;
    elseif size(EEG_raw.data(1,:),2) < srate*60*60 % less than an hour (in minutes)
        divideby = 60;
    else % more than an hour (in hours)
        divideby = 60*60;
    end
    % Number of channels
    num_channels = length(EEG_raw.chanlocs);

%     %%%%% Check the reference mastoid channels for their signal quality
%     % Step 1: Check for flatline signals (variance close to zero)
%     flatline_threshold = 1e-6;  % Threshold for detecting flatline channels
%     mastoid_variance = var(EEG_mastoids.data, 0, 2);  % Variance of each mastoid channel
%     flatline_channels = mastoid_variance < flatline_threshold;
%     if any(flatline_channels)
%         warning('One or more mastoid channels have a flatline signal.');
%     end
%     % Check their inter-mastoid-correlation
%     cor_mastoids = abs(corrcoef(EEG_mastoids.data(:,:)'));
%     correlation_value = cor_mastoids(1, 2); 
%     disp(['Correlation of Mastoid 1 with Mastoid 2: ', num2str(correlation_value)]);
%     % Check their individual correlation to the average of C3 and C4
%     % Find indices of C3 and C4 by their labels in EEG_raw
%     C3_index = find(strcmp({EEG_raw.chanlocs.labels}, 'C3'));
%     C4_index = find(strcmp({EEG_raw.chanlocs.labels}, 'C4'));
%     % Extract mastoid channels from EEG_mastoids (assuming they are channels 1 and 2)
%     mastoid1 = EEG_mastoids.data(1, :);
%     mastoid2 = EEG_raw.data(C3_index, :);
%     % Extract C3 and C4 channels using their labels
%     C3 = EEG_raw.data(C3_index, :);
%     C4 = EEG_raw.data(C4_index, :);
%     % Calculate the average of C3 and C4
%     avg_C3_C4 = mean([C3; C4]);
%     % Compute correlation between mastoid channels and the average of C3 and C4
%     corr_mastoid1 = corr(mastoid1', mastoid2');
%     corr_mastoid2 = corr(mastoid2', avg_C3_C4');
%     % Display the results
%     disp(['Correlation of Mastoid 1 with C3 : ', num2str(corr_mastoid1)]);
%     disp(['Correlation of Mastoid 2 with C3 & C4 average: ', num2str(corr_mastoid2)]);
% 
%     % Extract mastoid channels from EEG_mastoids (assuming they are channels 1 and 2)
%     mastoid1 = EEG_mastoids_nonoise.data(2, :);
%     mastoid2 = EEG_raw_nonoise.data(C3_index, :);
%     % Compute correlation between mastoid channels and the average of C3 and C4
%     corr_mastoid1 = corr(mastoid1', mastoid2');
%     disp(['Correlation of Mastoid 1 with C3 : ', num2str(corr_mastoid1)]);
% 
% 
%     % % Plot mastoids
%     % % Start of plotting the spectrogram
%     % fig = figure('visible', 'off', 'WindowStyle', 'normal', 'Color', 'w', 'OuterPosition', screenSize);
%     % % Maximize the figure window
%     % screenSize = get(0, 'ScreenSize'); % Get the screen size from the root window
%     % set(fig, 'Position', [1 1 screenSize(3) screenSize(4)]); % Set the figure size to cover the whole screen
%     % hold on;
%     % % Extract mastoid channel indices from EEG_mastoids
%     % channel_indices = 1:size(EEG_mastoids.data, 1);  % Assuming EEG_mastoids contains only the mastoid channels
%     %     for ith_channel_idx = 1:length(channel_indices)
%     %         ith_channel = channel_indices(ith_channel_idx);  % Use the specific channel index
%     %         ax = subplot(2, 1, ith_channel_idx);
%     %         spectrogram(EEG_mastoids.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
%     %         ylim(fpass);
%     %         caxis([-40, 40]);
%     %         colormap('jet');
%     %         cbar = colorbar;  % Create or get handle to the colorbar
%     %         cbar.Color = 'k';  % Set the color of the colorbar's text and ticks to black
%     %         cbar.TickDirection = 'out';  % Optionally set the ticks to point outwards
%     %         title(EEG_mastoids.chanlocs(ith_channel).labels, 'Interpreter', 'none');
%     % end
%     % % Add custom super title adjusted for better positioning
%     % annotation(fig, 'textbox', [0, 0.95, 1, 0.05], ...
%     %     'String', [study, ' ', participant, ' Mastoids, Correlation of ', num2str(correlation_value)], ...
%     %     'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
%     %     'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
%     % % Set text and axes colors to black and adjust title positions
%     % axs = findobj(fig, 'Type', 'axes');
%     % for ax = axs'
%     %     set(ax, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
%     %     set(ax.Title, 'Color', 'k');
%     % end
%     % txts = findobj(fig, 'Type', 'text');
%     % for txt = txts'
%     %     set(txt, 'Color', 'k');
%     % end
%     % hold off;
%     % % Saving the plot as PNG in pathtosave
%     % filename = [study, '_', participant, '_Mastoids.png'];
%     % print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
%     % % Display the full path of the saved file in the command window
%     % disp(['Saved ', filename, ' in ',  pathtosave]);
%     % close all hidden;
% 
% 
% 
% 
% 
% 
% % Combine mastoid channels from EEG_mastoids with C3 and C4 channels from EEG_raw
% channel_indices_mastoids = 1:size(EEG_mastoids.data, 1);  % Assuming EEG_mastoids has only the mastoid channels
% channel_indices_C3_C4 = [C3_index, C4_index];  % Indices of C3 and C4 in EEG_raw
% 
% % Plotting both mastoid and C3/C4 channels
% fig = figure('visible', 'off', 'WindowStyle', 'normal', 'Color', 'w', 'OuterPosition', screenSize);
% screenSize = get(0, 'ScreenSize'); % Get the screen size from the root window
% set(fig, 'Position', [1 1 screenSize(3) screenSize(4)]); % Set the figure size to cover the whole screen
% hold on;
% 
% % Plotting mastoid channels
% for ith_channel_idx = 1:length(channel_indices_mastoids)
%     ith_channel = channel_indices_mastoids(ith_channel_idx);  % Use the specific mastoid channel index
%     ax = subplot(4, 1, ith_channel_idx);  % Adjusted to 4 rows for both mastoids and C3/C4
%     spectrogram(EEG_mastoids.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
%     ylim(fpass);
%     caxis([-40, 40]);
%     colormap('jet');
%     cbar = colorbar;  % Create or get handle to the colorbar
%     cbar.Color = 'k';  % Set the color of the colorbar's text and ticks to black
%     cbar.TickDirection = 'out';  % Optionally set the ticks to point outwards
%     title(EEG_mastoids.chanlocs(ith_channel).labels, 'Interpreter', 'none');
% end
% 
% % Plotting C3 and C4 channels
% for ith_channel_idx = 1:length(channel_indices_C3_C4)
%     ith_channel = channel_indices_C3_C4(ith_channel_idx);  % Use the specific C3 or C4 channel index
%     ax = subplot(4, 1, ith_channel_idx + 2);  % Continue subplot from the next row (adjusting for 4 rows)
%     spectrogram(EEG_raw.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
%     ylim(fpass);
%     caxis([-40, 40]);
%     colormap('jet');
%     cbar = colorbar;  % Create or get handle to the colorbar
%     cbar.Color = 'k';  % Set the color of the colorbar's text and ticks to black
%     cbar.TickDirection = 'out';  % Optionally set the ticks to point outwards
%     title(EEG_raw.chanlocs(ith_channel).labels, 'Interpreter', 'none');
% end
% 
% % Add custom super title adjusted for better positioning
% annotation(fig, 'textbox', [0, 0.95, 1, 0.05], ...
%     'String', [study, ' ', participant, ' Mastoids CorrM1,M2 = ', num2str(correlation_value), ...
%     ', CorrM1,C3&C4 = ', num2str(corr_mastoid1), ', CorrM2,C3&C4 = ', num2str(corr_mastoid2)],...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
%     'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
% 
% % Set text and axes colors to black and adjust title positions
% axs = findobj(fig, 'Type', 'axes');
% for ax = axs'
%     set(ax, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
%     set(ax.Title, 'Color', 'k');
% end
% txts = findobj(fig, 'Type', 'text');
% for txt = txts'
%     set(txt, 'Color', 'k');
% end
% 
% hold off;
% 
% % Saving the plot as PNG in pathtosave
% filename = [study, '_', participant, '_Mastoids_C3_C4.png'];
% print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
% disp(['Saved ', filename, ' in ', pathtosave]);
% close all hidden;

























    %% bandpass and detrend, to remove sweat and other artifacts
    fprintf('---bandpass detrend\n');
    EEG_detrend = EEG_raw;
    if size(EEG_raw.data,1) > 0
        for chcha = 1:size(EEG_raw.data,1)
            EEG_detrend.data(chcha,:) = detrend(EEG_raw.data(chcha,:));
            EEG_detrend.data(chcha,:) = bandpass(EEG_detrend.data(chcha,:), [0.5, 40], EEG_raw.srate);
        end
    end
    

    %% detect bad signal to remove for bad channel detection
    fprintf('---detecting bad signal\n');
    [remove_channel_which, cut_chunks, rem_ind_all_matrix] = eeg_acrossfreq_artifact(EEG_detrend);
    rem_ind_all = sum(rem_ind_all_matrix,1) > 0; % von matrix zu vektor, falls irgendwo eine 1 ist (bad signal) wird der ganze Zeitpunkt auf 1 gesetzt
    

    %% remove electrical noise from sienna
    fprintf('---removal electrical noise\n');
    EEG_noise = EEG_detrend;
    EEG_noise = eeg_electricalnoise_reducer(EEG_detrend, rem_ind_all);


    %% detect bad channels without bad signal
    labels = {EEG_detrend.chanlocs.labels};
    labels = labels(ismember(labels, {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'}));
    bad_channels = labels(eeg_badchannel_detect(pop_select(EEG_noise, 'channel', {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'}), zeros(1,size(EEG_detrend.data,2)), 0.5))
    disp(['Bad channels: ', bad_channels]);


    %% interpolate bad channels
    % Check if bad_channels is empty
    if ~isempty(bad_channels)
        % Perform interpolation only if there are bad channels
        bad_channel_indices = eeg_chaninds(EEG_noise, bad_channels);
        EEG_interp = pop_interp(EEG_noise, bad_channel_indices);
    else
        % If bad_channels is empty, skip interpolation and keep EEG data unchanged
        EEG_interp = EEG_noise;
    end


    %% detect bad signal again to exclude
    fprintf('---detecting bad signal\n');
    [remove_channel_which, cut_chunks, rem_ind_all_matrix] = eeg_acrossfreq_artifact(EEG_interp);
    rem_ind_all = sum(rem_ind_all_matrix,1) > 0; % von matrix zu vektor, falls irgendwo eine 1 ist (bad signal) wird der ganze Zeitpunkt auf 1 gesetzt

   
    %% Calculate and plot summary of bad signal
    % Calculate the percentage of bad segments
    perc_bad = sum(rem_ind_all) / length(rem_ind_all) * 100;
    dur_bad = perc_bad / 100 * length(rem_ind_all) / EEG_detrend.srate / 60;
    % Create a new figure for the combined plot
    screenSize = get(0, 'ScreenSize');  % Returns a vector [left, bottom, width, height]
    fig = figure('visible', 'off', 'WindowStyle', 'normal', 'Color', 'w', 'OuterPosition', screenSize);
    % Calculate segment durations
    diff_vec = [0, diff(rem_ind_all), 0];  % Take the diff and pad with zeros at both ends
    segment_starts = find(diff_vec == 1);  % Start of each segment is where diff is 1
    segment_ends = find(diff_vec == -1) - 1;  % End of each segment is where diff is -1
    % Handle case where the last element is a 1 (extra start without end)
    if length(segment_starts) > length(segment_ends)
        % Add the last index as the end of the last segment if it ends with 1
        segment_ends = [segment_ends, length(rem_ind_all)];
    end
    % Handle case where the first element is a 1 
    if length(segment_starts) < length(segment_ends)
        % Add the first index as the beginnning 
        segment_starts = [1, segment_starts];
    end
    segment_lengths = segment_ends - segment_starts + 1;  % Calculate lengths of each segment in samples
    segment_durations = segment_lengths / EEG_interp.srate;  % Calculate duration of each segment in seconds
    % Left subplot: Histogram of segment durations
    subplot(1, 2, 1);
    histogram(segment_durations);
    title('Histogram of Segment Durations');
    xlabel('Duration (seconds)');
    ylabel('Count');
    % Calculate the total number of electrodes containing bad data at a given timepoint
    column_sums = sum(rem_ind_all_matrix);  % Sum across rows for each column (timepoint)
    nonzero_sums = column_sums(column_sums > 0);  % Filter out zero sums
    % Right subplot: Histogram of number of electrodes with bad data
    subplot(1, 2, 2);
    histogram(nonzero_sums, 'Normalization', 'count');
    title('Histogram of Bad Electrode Counts');
    xlabel('Number of Bad Electrodes');
    ylabel('Count');
    % Add a super title with the percentage and duration of bad segments
    % Add a big title (super title) with the overall summary
    big_title = sprintf([study, ' ', participant, ' Summary of bad data']);
    sgtitle(big_title, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');  % Big title
    % Add an additional annotation for the summary information
    annotation_text = sprintf('The percentage of bad segments is %.2f%%, which equals to %.2f minutes.', perc_bad, dur_bad);
    annotation('textbox', [0.5, 0.05, 0, 0], 'String', annotation_text, 'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'EdgeColor', 'none', 'Color', 'k');
    axs = findobj(fig, 'Type', 'axes');
    for ax = axs'
        set(ax, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
        set(ax.Title, 'Color', 'k');
    end
    txts = findobj(fig, 'Type', 'text');
    for txt = txts'
        set(txt, 'Color', 'k');
    end
    % Saving the plot as PNG in pathtosave
    filename = [study, '_', participant, '_Summary_of_bad_data.png'];
    print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
    % Display the full path of the saved file in the command window
    disp(['Saved ', filename, ' in ',  pathtosave]);
    close all hidden;

    % Generate and save the bad data
    EEG_bad = EEG_interp;
    EEG_bad.data = [];
    EEG_bad.times = [];
    EEG_bad.data = EEG_interp.data(:, rem_ind_all);
    EEG_bad.times = EEG_interp.times(:, rem_ind_all);
    filePath = fullfile(pathtosave, [study, '_', participant, '_EEG_bad_data.mat']);
    % Save the EEG_EOG variable in the specified path
    save(filePath, 'EEG_bad');
    % Optionally display the path where the file was saved
    disp(['Saved bad data in ', filePath]);


    %% remove bad signal from signal
    EEG_interp_removed = EEG_interp;
    % Create a logical index for good data (where rem_ind_all is 0)
    good_data_indices = ~rem_ind_all;  % The tilde (~) operator negates the logical array
    % Filter the entire EEG_detrend_removed.data and .times matrix to remove bad data points
    % This keeps only the columns (time points) where rem_ind_all is 0 (or false for good_data_indices)
    EEG_interp_removed.data = EEG_interp.data(:, good_data_indices);
    EEG_interp_removed.times = EEG_interp.times(:, good_data_indices);
    ECG.data = ECG.data(:,good_data_indices);
    ECG.times = ECG.times(:,good_data_indices);

    %% Rereference to AR (average reference)
    EEG_reref = pop_reref(EEG_raw, []);
    %% Rereference using the REST algorithm
    EEG = EEG_raw;
    eeglab redraw;

    

    %% remove heart beat artifacts here again I require the movement artifacts to be excluded (eeg_jelena_ecg_removal will probably not yet run correctly, I will need to adjust it)
    fprintf('---ECG removal\n');
    EEG_ECG = EEG_noise;
    EEG_ECG = eeg_jelena_ecg_removal(EEG_noise, ECG, cut_chunks);


    %% let ASR run sequentially in windows
    fprintf('---ASR\n');
        overlapDuration = windowDuration / 2;
        windowSamples = windowDuration * 60 * srate;  % convert window duration to samples
        overlapSamples = overlapDuration * 60 * srate;  % convert overlap duration to samples
        nSamples = size(EEG_ECG.data, 2);  % total number of samples in the EEG data
        segmentStarts = 1:overlapSamples:(nSamples - windowSamples + 1);
        segmentEnds = segmentStarts + windowSamples - 1;
        % Adjust the last segment to end precisely at the end of the data
        if segmentEnds(end) > nSamples
            segmentEnds(end) = nSamples;
        end
        % Check if the last segment is less than the standard window duration
        if (segmentEnds(end) - segmentStarts(end) + 1) < windowSamples
            % If the last segment is too short, extend the previous segment
            % Remove the last start and set the second-to-last end to the end of data
            segmentStarts(end) = [];
            segmentEnds(end-1) = nSamples;
        end
        % Adjust last segment's end again after modification to cover the end of data
        if segmentEnds(end) < nSamples
            segmentEnds(end) = nSamples;
        end
        % delete possibly existing jobs in parallel pool
        myCluster = parcluster('Processes');
        delete(myCluster.Jobs);
        % Open a parallel pool if it isn't already open
        if isempty(gcp('nocreate'))
            parpool;  % This opens a parallel pool with the default number of workers
        end
        % Sequentially apply ASR (this takes a while)
        EEG_segments = cell(1, numel(segmentStarts));
        parfor j = 1:numel(segmentStarts)
            segmentData = EEG_ECG.data(:, segmentStarts(j):segmentEnds(j));
            % Create a temporary EEG structure for ASR
            EEG_temp = EEG_ECG;
            EEG_temp.data = segmentData;
            EEG_temp.pnts = size(segmentData, 2);
            EEG_temp.xmax = EEG_temp.xmin + (EEG_temp.pnts - 1) / EEG_temp.srate;
            % Apply ASR
            [~, ~, EEG_clean, ~] = clean_artifacts(EEG_temp, ...
                'Highpass', 'off', ...
                'ChannelCriterion', 'off', ...
                'FlatlineCriterion', 'off', ...
                'LineNoiseCriterion', 'off', ...
                'BurstCriterion', SDCutoff, ...
                'WindowCriterion', 0.25);
            EEG_segments{j} = EEG_clean.data;
        end
        
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end

        % Create the tapers to smoothen the overlaps of the data segments
        EEG_segments_tapered = cell(1, numel(segmentStarts));
        for j = 1:numel(EEG_segments)
            segmentLength = size(EEG_segments{j}, 2);  % Current segment length
            taperLength = overlapSamples;  % Length of the taper
            if j == 1  % First segment
                % Starts full, tapers down at the end
                taper = [ones(1, segmentLength - taperLength), cos(linspace(0, pi/2, taperLength)).^2];
            elseif j == numel(EEG_segments)  % Last segment
                % Tapers up at the start, remains full for the rest
                taper = [cos(linspace(pi/2, 0, taperLength)).^2, ones(1, segmentLength - taperLength)];
            else  % Middle segments
                % Tapers up at the start, full in the middle, tapers down at the end
                taper = [cos(linspace(pi/2, 0, taperLength)).^2, ones(1, segmentLength - 2 * taperLength), cos(linspace(0, pi/2, taperLength)).^2];
            end
            % Apply the taper
            EEG_segments_tapered{j} = EEG_segments{j} .* repmat(taper, size(EEG_segments{j}, 1), 1);
        end
        
        % Put the tapered segments back together
        EEG_EOG = EEG_ECG;
        EEG_EOG.data = zeros(size(EEG_EOG.data));  % Initialize with zeros for addition
        currentSample = 1;
        for j = 1:numel(EEG_segments_tapered)
            segmentLength = size(EEG_segments_tapered{j}, 2);  % Determine the length of the current segment
            % Calculate the end index for this segment
            endIndex = currentSample + segmentLength - 1;
            % Add the segment to the existing data in EEG_EOG, summing where segments overlap
            if endIndex > size(EEG_EOG.data, 2)
                endIndex = size(EEG_EOG.data, 2);  % Ensure not to exceed the bounds
                segmentLength = endIndex - currentSample + 1;  % Adjust segment length accordingly
            end
            EEG_EOG.data(:, currentSample:endIndex) = ...
                EEG_EOG.data(:, currentSample:endIndex) + EEG_segments_tapered{j}(:, 1:segmentLength);
            % Update the currentSample index for the next segment, shifting by the overlap
            currentSample = currentSample + overlapSamples;  % Move by the overlap amount
        end

    %% Calculate spectrum for each step
    % For EEG_raw (Raw data)
    fprintf('---calculate the spectrum of every preprocessing step\n');
    [spectra1, freqs1] = spectopo(EEG_raw.data, 0, EEG_raw.srate, 'plot', 'off');
    % For EEG_detrend
    [spectra2, freqs2] = spectopo(EEG_detrend.data, 0, EEG_detrend.srate, 'plot', 'off');
    % For EEG_noise
    [spectra3, freqs3] = spectopo(EEG_noise.data, 0, EEG_noise.srate, 'plot', 'off');
    % For EEG_ECG
    [spectra4, freqs4] = spectopo(EEG_ECG.data, 0, EEG_ECG.srate, 'plot', 'off');
    % For EEG_EOG
    [spectra5, freqs5] = spectopo(EEG_EOG.data, 0, EEG_EOG.srate, 'plot', 'off');
   
  


%% Start of plotting the spectrogram
fig = figure('visible', 'off', 'WindowStyle', 'normal', 'Color', 'w', 'OuterPosition', screenSize);
% Maximize the figure window
screenSize = get(0, 'ScreenSize'); % Get the screen size from the root window
set(fig, 'Position', [1 1 screenSize(3) screenSize(4)]); % Set the figure size to cover the whole screen
hold on;
% Array of EEG structures and labels
EEG_steps = {EEG_raw, EEG_detrend, EEG_interp_removed, EEG_noise, EEG_ECG, EEG_EOG};
step_labels = {'Raw data', 'Bandpass & Detrend & bad signal', 'Bad signal removal', 'Electrical noise removal', 'ECG removal', 'ASR'};
% Define the array of cut chunks for each step if applicable
cut_chunks_steps = {[], cut_chunks, [], [], [], []};  % Adjust according to when chunks are used
% List of specific channels to display
specific_channels = {'F3', 'C3', 'P3', 'O1'};
% Find indices of these channels in the EEG data
channel_indices = find(ismember({EEG_raw.chanlocs.labels}, specific_channels));
% Iterate through each preprocessing step
for step = 1:length(EEG_steps)
    for ith_channel_idx = 1:length(channel_indices)
        ith_channel = channel_indices(ith_channel_idx);  % Use the specific channel index
        ax = subplot(length(EEG_steps), length(channel_indices), ith_channel_idx + (step - 1) * length(channel_indices));
        spectrogram(EEG_steps{step}.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
        ylim(fpass);
        caxis([-40, 40]);
        colormap('jet');
        cbar = colorbar;  % Create or get handle to the colorbar
        cbar.Color = 'k';  % Set the color of the colorbar's text and ticks to black
        cbar.TickDirection = 'out';  % Optionally set the ticks to point outwards
        title(EEG_steps{step}.chanlocs(ith_channel).labels, 'Interpreter', 'none');
        % Add black bars if cut chunks are available for this step
        if ~isempty(cut_chunks_steps{step}) && ~isempty(cut_chunks_steps{step}{ith_channel})
            ccunk = cut_chunks_steps{step}{ith_channel};
            for ccc = 1:size(ccunk, 1)
                x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                patch(x, y, 'black', 'FaceAlpha', 1);
            end
        end
    end
    % Position of the current axes
    pos = get(ax, 'Position');
    % Create an annotation that spans the full width of the figure
    annotation('textbox', [0, pos(2) + pos(4), 1, 0.05], ...
               'String', step_labels{step}, ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'Color', 'k');
end
% Add custom super title adjusted for better positioning
annotation(fig, 'textbox', [0, 0.95, 1, 0.05], ...
           'String', [study, ' ', participant, ' spectrogram'], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
           'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
% Set text and axes colors to black and adjust title positions
axs = findobj(fig, 'Type', 'axes');
for ax = axs'
    set(ax, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    set(ax.Title, 'Color', 'k');
end
txts = findobj(fig, 'Type', 'text');
for txt = txts'
    set(txt, 'Color', 'k');
end
hold off;
% Saving the plot as PNG in pathtosave
filename = [study, '_', participant, '_Spectrogram.png'];
print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
% Display the full path of the saved file in the command window
disp(['Saved ', filename, ' in ',  pathtosave]);
close all hidden;







% Create a figure for all power spectra
fprintf('---plot the individual spectra\n');
fig = figure('visible', 'off', 'Color', 'w', 'OuterPosition', screenSize);
% Plot the spectra
plot(freqs1, spectra1(1,:), 'Marker', 'o', 'LineStyle', '--', 'Color', '#003049');
hold on;
%plot(freqs2, spectra2(1,:), 'Marker', '*', 'LineStyle', '--', 'Color', '#369279');
%plot(freqs3, spectra3(1,:), 'Marker', '.', 'LineStyle', '--', 'Color', '#55CC95');
%plot(freqs4, spectra4(1,:), 'Marker', '<', 'LineStyle', '--', 'Color', '#2B8B5E');
plot(freqs5, spectra5(1,:), 'Marker', '*', 'LineStyle', '--', 'Color', '#2B8B5E');
ax = gca;  % Gets handle to the current axes
set(ax, 'XColor', 'k', 'YColor', 'k');  % Sets the x-axis and y-axis color to black
set(ax, 'Color', 'w');  % Optional: Explicitly set axes background color if needed
% Setting labels, title, and legend
xlabel(ax, 'Frequency (Hz)', 'Color', 'k');  % Ensure text is black
ylabel(ax, 'Power (dB)', 'Color', 'k');
title(ax, 'Power Spectra of raw and preprocessed data', 'Color', 'k');
legend(ax, {'Raw', 'Preprocessed'}, 'TextColor', 'k', 'Color', 'w');
%legend(ax, {'Raw data', 'Detrend', 'Noise', 'ECG', 'ASR'}, 'TextColor', 'k', 'Color', 'w');
% Set grid and limits
grid on;
xlim(ax, [0 40]);
% Save the plot
filename = [study, '_', participant, '_Powerspectrum_individual.png'];
print(fig, fullfile(pathtosave, filename), '-dpng', '-r300');
disp(['Saved ', filename, ' in ', pathtosave]);
% Release the hold on the plot
hold off;
close all hidden;





% Create a figure for all spectrogram-differences
fprintf('---plot the differences of the spectra\n');
figure('visible', 'off', 'Color', 'w', 'OuterPosition', screenSize);
% Plot the first spectrum
plot(freqs1, spectra1(1,:)-spectra2(1,:), 'Marker', 'o', 'LineStyle', '--', 'Color', '#003049');
hold on;
plot(freqs2, spectra2(1,:)-spectra3(1,:), 'Marker', '.', 'LineStyle', '--', 'Color', '#16575C');
hold on;
plot(freqs3, spectra3(1,:)-spectra4(1,:), 'Marker', '*', 'LineStyle', '--', 'Color', '#55CC95');
hold on;
plot(freqs3, spectra4(1,:)-spectra5(1,:), 'Marker', 'square', 'LineStyle', '--', 'Color', '#004A26');
% Adding labels and title
ax = gca;  % Gets handle to the current axes
set(ax, 'XColor', 'k', 'YColor', 'k');  % Sets the x-axis and y-axis color to black
set(ax, 'Color', 'w');  % Optional: Explicitly set axes background color if needed
% Setting labels, title, and legend
xlabel(ax, 'Frequency (Hz)', 'Color', 'k');  % Ensure text is black
ylabel(ax, 'Power (dB)', 'Color', 'k');
title(ax, 'Power Spectra difference of preprocessing steps', 'Color', 'k');
legend(ax, {'Raw - Detrend', 'Detrend - Noise', 'Noise - ECG', 'ECG - ASR'}, 'TextColor', 'k', 'Color', 'w', 'Location', 'northwest');
% Enable grid for better visibility
grid on;
xlim(ax, [0 40]);
% Release the hold on the plot
hold off;
% Saving the plot as png in pathtosave
filename = [study, '_', participant, '_Powerspectrum_difference.png'];
print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
% Display the full path of the saved file in the command window
disp(['Saved ', filename, ' in ',  pathtosave]);
close all hidden;




filePath = fullfile(pathtosave, [study, '_', participant, '_EEG_preprocessed.mat']);
% Save the EEG_EOG variable in the specified path
save(filePath, 'EEG_EOG');
% Optionally display the path where the file was saved
disp(['Saved EEG_EOG in ', filePath]);

clearvars -except folderpath i edf_files web_files study selected_channels windowDuration SDCutoff
end