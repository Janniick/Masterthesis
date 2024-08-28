% other functions (maybe here not needed, but for feature extraction at least they are needed
addpath(genpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_mesmart\scripts'));

% preprocessing Funktionen
addpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\matlab_processing_benji\');

% eeglab path
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui; % Initialize EEGLAB (eeglab nogui; if you want no GUI)

% plot Funktion
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing');

% edf & web file directory
folderpath = 'D:\Masterarbeit Jannick\Data\GHB_TRA\GHB_TRA_EEG';
edf_files = dir(fullfile(folderpath, '*\*.edf'));
web_files = dir(fullfile(folderpath, '*\*.web'));

% loop through all files
for i = 1:numel(edf_files) %randperm(numel(edf_files)) % 1:numel(edf_files)
    fprintf('********************************************%s****************************************\n', repmat('*', 1, length(edf_files(i).name)));
    fprintf('---------------------------------------File:%s----------------------------------------\n', edf_files(i).name);
    fprintf('                                  Iteration: %s\n', string(i));
    fprintf('                             Percentage run: %s\n', strcat(string(round(100*i/numel(edf_files),0)), '%'));
    fprintf('********************************************%s****************************************\n', repmat('*', 1, length(edf_files(i).name)));
    
    % which filepath
    edf_name = edf_files(i).name;
    edf_path = fullfile(edf_files(i).folder, edf_name);
    web_name = web_files(i).name;
    web_path = fullfile(web_files(i).folder, web_name);
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
    ECG = pop_select(EEG_raw, 'channel', {'ECG'});
    EEG_raw = pop_select(EEG_raw, 'channel', {'E1', 'E2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', });
    %EEG = pop_select(EEG_raw, 'channel', {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'});
    newchannelLabels = {EEG_raw.chanlocs.labels};

    % downsample
    fprintf('---downsample\n');
    EEG_raw = pop_resample(EEG_raw, 256);
    ECG = pop_resample(ECG, 256);
    
    %%%%%%%%%%% Define start and end times in minutes for the preffered window
    condition = 'fullnightremoved';

    % define parameters for plotting
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


    %%%%% bandpass and detrend, to remove sweat and other artifacts
    fprintf('---bandpass detrend\n');
    EEG_detrend = EEG_raw;
    if size(EEG_raw.data,1) > 0
        for chcha = 1:size(EEG_raw.data,1)
            EEG_detrend.data(chcha,:) = detrend(EEG_raw.data(chcha,:));
            EEG_detrend.data(chcha,:) = bandpass(EEG_detrend.data(chcha,:), [0.5, 40], EEG_raw.srate);
        end
    end
    

    %%%%% detect bad signal to remove from signal
    fprintf('---detecting bad signal\n');
    [remove_channel_which, cut_chunks, rem_ind_all] = eeg_acrossfreq_artifact(EEG_detrend);
    rem_ind_all = sum(rem_ind_all,1) > 0; % von matrix zu vektor, falls irgendwo eine 1 ist (bad signal) wird der ganze Zeitpunkt auf 1 gesetzt
    

    %%%%% remove bad signal from signal
    EEG_detrend_removed = EEG_detrend;
    % Create a logical index for good data (where rem_ind_all is 0)
    good_data_indices = ~rem_ind_all;  % The tilde (~) operator negates the logical array
    % Filter the entire EEG_detrend_removed.data and .times matrix to remove bad data points
    % This keeps only the columns (time points) where rem_ind_all is 0 (or false for good_data_indices)
    EEG_detrend_removed.data = EEG_detrend.data(:, good_data_indices);
    EEG_detrend_removed.times = EEG_detrend.times(:, good_data_indices);
    ECG.data = ECG.data(:,good_data_indices);
    ECG.times = ECG.times(:,good_data_indices);

   
    %%%%% remove electrical noise from sienna, and make sure to exclude the movement artifacts for calculation (output EEG still has all data, including the movments)
    fprintf('---removal electrical noise\n');
    EEG_noise = EEG_detrend;
    EEG_noise = eeg_electricalnoise_reducer(EEG_detrend, rem_ind_all);
    EEG_noise.data = EEG_noise.data(:, good_data_indices);
    EEG_noise.times = EEG_noise.times(:, good_data_indices);
    

    %%%%% remove heart beat artifacts here again I require the movement artifacts to be excluded (eeg_jelena_ecg_removal will probably not yet run correctly, I will need to adjust it)
    fprintf('---ECG removal\n');
    EEG_ECG = EEG_noise;
    EEG_ECG = eeg_jelena_ecg_removal(EEG_noise, ECG, cut_chunks);


    %%%%% let the ASR run sequentially in windows of x minutes
    fprintf('---ASR\n');
        windowDuration = 8; % window length in minutes
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
        % Open a parallel pool if it isn't already open
        if isempty(gcp('nocreate'))
            parpool;  % This opens a parallel pool with the default number of workers
        end
        % Sequentially apply ASR (this takes a while)
        EEG_segments = cell(1, numel(segmentStarts));
        parfor i = 1:numel(segmentStarts)
            segmentData = EEG_ECG.data(:, segmentStarts(i):segmentEnds(i));
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
                'BurstCriterion', 30, ...
                'WindowCriterion', 0.25);
            EEG_segments{i} = EEG_clean.data;
        end
        
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end

        % Create the tapers to smoothen the overlaps of the data segments
        EEG_segments_tapered = cell(1, numel(segmentStarts));
        for i = 1:numel(EEG_segments)
            segmentLength = size(EEG_segments{i}, 2);  % Current segment length
            taperLength = overlapSamples;  % Length of the taper
            if i == 1  % First segment
                % Starts full, tapers down at the end
                taper = [ones(1, segmentLength - taperLength), cos(linspace(0, pi/2, taperLength)).^2];
            elseif i == numel(EEG_segments)  % Last segment
                % Tapers up at the start, remains full for the rest
                taper = [cos(linspace(pi/2, 0, taperLength)).^2, ones(1, segmentLength - taperLength)];
            else  % Middle segments
                % Tapers up at the start, full in the middle, tapers down at the end
                taper = [cos(linspace(pi/2, 0, taperLength)).^2, ones(1, segmentLength - 2 * taperLength), cos(linspace(0, pi/2, taperLength)).^2];
            end
            % Apply the taper
            EEG_segments_tapered{i} = EEG_segments{i} .* repmat(taper, size(EEG_segments{i}, 1), 1);
        end
        
        % Put the tapered segments back together
        EEG_EOG = EEG_ECG;
        EEG_EOG.data = zeros(size(EEG_EOG.data));  % Initialize with zeros for addition
        currentSample = 1;
        for i = 1:numel(EEG_segments_tapered)
            segmentLength = size(EEG_segments_tapered{i}, 2);  % Determine the length of the current segment
            % Calculate the end index for this segment
            endIndex = currentSample + segmentLength - 1;
            % Add the segment to the existing data in EEG_EOG, summing where segments overlap
            if endIndex > size(EEG_EOG.data, 2)
                endIndex = size(EEG_EOG.data, 2);  % Ensure not to exceed the bounds
                segmentLength = endIndex - currentSample + 1;  % Adjust segment length accordingly
            end
            EEG_EOG.data(:, currentSample:endIndex) = ...
                EEG_EOG.data(:, currentSample:endIndex) + EEG_segments_tapered{i}(:, 1:segmentLength);
            % Update the currentSample index for the next segment, shifting by the overlap
            currentSample = currentSample + overlapSamples;  % Move by the overlap amount
        end

    % Calculate spectrum for each step
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
   
  


% Start of plotting the spectrogram
fig = figure('visible', 'on', 'WindowStyle', 'normal', 'Color', 'w');
% Maximize the figure window
screenSize = get(0, 'ScreenSize'); % Get the screen size from the root window
set(fig, 'Position', [1 1 screenSize(3) screenSize(4)]); % Set the figure size to cover the whole screen
hold on;
% Array of EEG structures and labels
EEG_steps = {EEG_raw, EEG_detrend, EEG_detrend_removed, EEG_noise, EEG_ECG, EEG_EOG};
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
% sgtitle([participant, ' ', condition, ' spectrogram'], 'FontSize', 14);
% Add custom super title adjusted for better positioning
annotation(fig, 'textbox', [0, 0.95, 1, 0.05], ...
           'String', [participant, ' ', condition, ' spectrogram'], ...
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
filename = [participant, '_', condition, '_Spectrogram.png'];
print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
% Display the full path of the saved file in the command window
disp(['Saved ', filename, ' in ',  pathtosave]);






% Create a figure for all power spectra
fprintf('---plot the individual spectra\n');
fig = figure('visible', 'off', 'Color', 'w');
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
filename = [participant, '_', condition, '_Powerspectrum_individual.png'];
print(fig, fullfile(pathtosave, filename), '-dpng', '-r300');
disp(['Saved ', filename, ' in ', pathtosave]);
% Release the hold on the plot
hold off;





% Create a figure for all spectrogram-differences
fprintf('---plot the differences of the spectra\n');
figure('visible', 'off', 'Color', 'w');
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
filename = [participant, '_', condition, '_Powerspectrum_difference.png'];
print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
% Display the full path of the saved file in the command window
disp(['Saved ', filename, ' in ',  pathtosave]);




filePath = fullfile(pathtosave, [participant, '_EEG_preprocessed.mat']);
% Save the EEG_EOG variable in the specified path
save(filePath, 'EEG_EOG');
% Optionally display the path where the file was saved
disp(['Saved EEG_EOG in ', filePath]);
end








    % Initialize a new EEG structure
% selectedEEG = EEG_ECG;
% 
% % Identify the indices of the channels you need
% pre_channels = {'E1', 'E2', 'F3', 'F4'};
% diff_channels = {'F3', 'F4'};
% pre_idx = find(ismember({EEG_ECG.chanlocs.labels}, pre_channels));
% diff_idx = find(ismember({EEG_ECG.chanlocs.labels}, diff_channels));
% 
% % Calculate the difference data for the specified channels
% difference_data = EEG_EOG.data(diff_idx, :) - EEG_ECG.data(diff_idx, :);
% 
% % Combine the selected pre-processed data with the difference data
% selectedEEG.data = [EEG_ECG.data(pre_idx, :); difference_data];
% 
% % Correct number of channels to include
% num_pre_channels = length(pre_idx);
% num_diff_channels = length(diff_idx);  % This will be less or equal to num_pre_channels depending on overlap
% 
% % Initialize a new chanlocs array with the correct number of entries
% new_chanlocs(num_pre_channels + num_diff_channels) = EEG_ECG.chanlocs(1);  % Copy the structure template
% 
% % Fill in the new chanlocs array
% idx = 1;
% for i = 1:num_pre_channels
%     new_chanlocs(idx) = EEG_ECG.chanlocs(pre_idx(i));
%     idx = idx + 1;
% end
% for i = 1:num_diff_channels
%     new_chanlocs(idx) = EEG_ECG.chanlocs(diff_idx(i));
%     new_chanlocs(idx).labels = [new_chanlocs(idx).labels '_diff'];  % Append '_diff' to the label
%     idx = idx + 1;
% end
% 
% % Assign the newly created chanlocs to selectedEEG
% selectedEEG.chanlocs = new_chanlocs;
% 
% % Update the number of channels in selectedEEG
% selectedEEG.nbchan = length(new_chanlocs);  % Should now reflect the correct count
% 
% % Visualize the new EEG structure
% pop_eegplot(selectedEEG, 1, 1, 0);