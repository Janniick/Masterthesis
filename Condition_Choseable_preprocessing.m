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

% omit bad files (if necessary)
% cut_omits = ismember({edf_files.name}, {'MGT_H17_N4_sleep_sienna2_732r0s1.edf', 'MGT_H23_N4_sleep_sienna2_808r0s2.edf', 'MGT_H23_N4_sleep_sienna2_808r0s3.edf'});
% edf_files = edf_files(setdiff(1:size(edf_files,1), find(cut_omits)));


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
    % Read the corresponding scoring .web file
    Scoring = fileread(web_path);

    % select channels
    channelLabels = {EEG_raw.chanlocs.labels};
    ECG = pop_select(EEG_raw, 'channel', {'ECG'});
    EEG_raw = pop_select(EEG_raw, 'channel', {'F3', 'F4'});
    %EEG = pop_select(EEG_raw, 'channel', {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'});
    newchannelLabels = {EEG_raw.chanlocs.labels};

    % downsample
    fprintf('---downsample\n');
    EEG_raw = pop_resample(EEG_raw, 256);
    ECG = pop_resample(ECG, 256);
    
    % Define start and end times in minutes for the preffered window
    condition = 'fullnight';
    % startTime = 2 * 60 + 0;  % Start at 1 hours and 0 minutes
    % endTime = startTime + 16;  % End time 16 minutes later
    % % Convert times to seconds
    % startSecond = startTime * 60;
    % endSecond = endTime * 60;
    % % Calculate sample indices
    % startIndex = round(startSecond * EEG_raw.srate) + 1;
    % endIndex = round(endSecond * EEG_raw.srate);
    % % Extract EEG data for the chosen window
    % selectedEEGData = EEG_raw.data(:, startIndex:endIndex);
    % % Create a new EEG structure for the chosen window
    % EEG_raw.data = selectedEEGData;  % Replace data with the chosen window
    % EEG_raw.pnts = size(selectedEEGData, 2);  % Update the number of points
    % EEG_raw.xmax = EEG_raw.xmin + (EEG_raw.pnts - 1) / EEG_raw.srate;  % Update xmax accordingly
    % % Extract ECG data for the chosen window
    % selectedECGData = ECG.data(:, startIndex:endIndex);
    % % Create a new ECG structure for the selected window
    % ECG.data = selectedECGData;  % Replace data with the selected window
    % ECG.pnts = size(selectedECGData, 2);  % Update the number of points
    % ECG.xmax = ECG.xmin + (ECG.pnts - 1) / ECG.srate;  % Update xmax accordingly

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
    % Number of channels and preprocessing steps
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
    
   
    %%%%% remove electrical noise from sienna, and make sure to exclude the movement artifacts for calculation (output EEG still has all data, including the movments)
    fprintf('---removal electrical noise\n');
    EEG_noise = EEG_detrend;
    rem_ind_all = sum(rem_ind_all,1) > 0; % von matrix zu vektor, falls irgendwo eine 1 ist (bad signal) wird der ganze Zeitpunkt auf 1 gesetzt
    EEG_noise = eeg_electricalnoise_reducer(EEG_detrend, rem_ind_all);
    

    %%%%% remove heart beat artifacts here again I require the movement artifacts to be excluded (eeg_jelena_ecg_removal will probably not yet run correctly, I will need to adjust it)
    fprintf('---ECG removal\n');
    EEG_ECG = EEG_noise;
    EEG_ECG = eeg_jelena_ecg_removal(EEG_noise, ECG, cut_chunks);


    %%%%% do some EOG removal here
    fprintf('---ASR\n');
    EEG_EOG = EEG_ECG;
    [~,~,EEG_EOG,removed_channels] = clean_artifacts(EEG_EOG, ...
        'Highpass','off', ...
        'ChannelCriterion', 'off', ... % 0 or 0.3 correlation? high density EEG is different
        'FlatlineCriterion', 'off', ...
        'LineNoiseCriterion', 'off', ...
        'BurstCriterion', 30, ... %burst criterion!!
        'WindowCriterion', 0.25);
    warning('off', 'MATLAB:colon:nonIntegerIndex');


    % Calculate spectrum for each step
    % For EEG_REM (Raw data)
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
fig = figure('visible', 'off', 'WindowStyle', 'normal');
% Maximize the figure window
screenSize = get(0, 'ScreenSize'); % Get the screen size from the root window
set(fig, 'Position', [1 1 screenSize(3) screenSize(4)]); % Set the figure size to cover the whole screen
hold on;
% Array of EEG structures and labels
EEG_steps = {EEG_raw, EEG_detrend, EEG_noise, EEG_ECG, EEG_EOG, EEG_EOG};
step_labels = {['Raw data of ', condition], 'Bandpass & Detrend', 'Electrical noise removal', 'ECG removal', 'EOG removal', 'Bad signal removal'};
% Define the array of cut chunks for each step if applicable
cut_chunks_steps = {[], [], [], [], [], cut_chunks};  % Adjust according to when chunks are used
% Iterate through each preprocessing step
for step = 1:length(EEG_steps)
    for ith_channel = 1:num_channels
        % Define subplot position and get handle
        ax = subplot(length(EEG_steps), num_channels, ith_channel + (step - 1) * num_channels);
        % Plot spectrogram for current channel and step
        spectrogram(EEG_steps{step}.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
        ylim(fpass);
        caxis([-40, 40]);
        colormap('jet');
        colorbar;
        % Title for each channel
        title(EEG_raw.chanlocs(ith_channel).labels, 'Interpreter', 'none');
        % Add black bars if cut chunks are available for this step
        if ~isempty(cut_chunks_steps{step})
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
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12, 'VerticalAlignment', 'bottom');
end
% Release hold after plotting all subplots
hold off;
% Saving the plot as PNG in pathtosave
filename = [participant, '_', condition, '_Spectrogram.png'];
print(gcf, fullfile(pathtosave, filename), '-dpng', '-r300');
% Display the full path of the saved file in the command window
disp(['Saved ', filename, ' in ',  pathtosave]);



 % Create a figure for all power spectra
    fprintf('---plot the individual spectra\n');
    figure; %('visible', 'on');
    % Plot the first spectrum 
    plot(freqs1, spectra1(1,:), 'blue-o');
    hold on; 
    plot(freqs2, spectra2(1,:), 'cyan-*'); 
    hold on; 
    plot(freqs3, spectra3(1,:), 'magenta-o'); 
    hold on; 
    plot(freqs4, spectra4(1,:), 'red-*'); 
    hold on; 
    plot(freqs5, spectra5(1,:), 'green-o'); 
    % Adding labels and title
    xlabel('Frequency (Hz)');
    xlim([0 40]);
    ylabel('Power (dB)');
    title(['Individual Power Spectra during ', condition]);
    % Add a legend to identify which line corresponds to which EEG data
    % legend({'Raw data', 'Detrend', 'Noise', 'ECG', 'EOG'});
    legend({'M1', 'M2', 'A1', 'A2', 'C3'});
    % Enable grid for better visibility
    grid on;
    % Release the hold on the plot
    hold off;
    % Saving the plot as png in pathtosave
    filename = [participant, '_', condition, '_Powerspectrum_individual.png'];
    print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
    % Display the full path of the saved file in the command window
    disp(['Saved ', filename, ' in ',  pathtosave]);


     % Create a figure for all spectrogram-differences
    fprintf('---plot the differences of the spectra\n');
    figure('visible', 'off');
    % Plot the first spectrum 
    plot(freqs1, spectra1(1,:)-spectra2(1,:), 'blue-o');
    hold on; 
    plot(freqs2, spectra2(1,:)-spectra3(1,:), 'cyan-*'); 
    hold on; 
    plot(freqs3, spectra3(1,:)-spectra4(1,:), 'magenta-o'); 
    hold on; 
    plot(freqs3, spectra4(1,:)-spectra5(1,:), 'green-*'); 
    % Adding labels and title
    xlabel('Frequency (Hz)');
    xlim([0 40]);
    ylabel('Power (dB)');
    title(['Comparison of EEG Power Spectra during ', condition]);
    % Add a legend to identify which line corresponds to which EEG data
    legend({'Raw data - Detrend', 'Detrend - Noise', 'Noise - ECG', 'ECG - EOG'});
    % Enable grid for better visibility
    grid on;
    % Release the hold on the plot
    hold off;
    % Saving the plot as png in pathtosave
    filename = [participant, '_', condition, '_Powerspectrum_difference.png'];
    print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
    % Display the full path of the saved file in the command window
    disp(['Saved ', filename, ' in ',  pathtosave]);

    savename = fullfile(pathtosave, [name, '_preprocessed.mat']);
    save(savename, 'EEG_EOG');% later in another feature extraction script you could load('EEG_preprocessed.mat'), of course by including the edf_name in the name
    % Display the full path of the saved file in the command window
    disp(['Saved ', savename, ' in ',  pathtosave]);
end








%%%%%%%%%%%%% Check mastoids

% Find the index of 'M1' in the labels array
indexM1 = find(strcmp({EEG_raw.chanlocs.labels},'M1'))
indexM2 = find(strcmp({EEG_raw.chanlocs.labels},'M2'))
indexA1 = find(strcmp({EEG_raw.chanlocs.labels},'A1'))
indexA2 = find(strcmp({EEG_raw.chanlocs.labels},'A2'))
indexC3 = find(strcmp({EEG_raw.chanlocs.labels},'C3'))

    [spectra1, freqs1] = spectopo(EEG_raw.data(indexM1,:), 0, EEG_raw.srate, 'plot', 'off');
    % For EEG_detrend
    [spectra2, freqs2] = spectopo(EEG_raw.data(indexM2,:), 0, EEG_raw.srate, 'plot', 'off');
    % For EEG_noise
    [spectra3, freqs3] = spectopo(EEG_raw.data(indexA1,:), 0, EEG_raw.srate, 'plot', 'off');
    % For EEG_ECG
    [spectra4, freqs4] = spectopo(EEG_raw.data(indexA2,:), 0, EEG_raw.srate, 'plot', 'off');
    % For EEG_EOG
    [spectra5, freqs5] = spectopo(EEG_raw.data(indexC3,:), 0, EEG_raw.srate, 'plot', 'off');

    fprintf('---plot the individual spectra\n');
    figure; %('visible', 'on');
    % Plot the first spectrum 
    plot(freqs1, spectra1(1,:), 'blue-o');
    hold on; 
    plot(freqs2, spectra2(1,:), 'cyan-*'); 
    hold on; 
    plot(freqs3, spectra3(1,:), 'magenta-o'); 
    hold on; 
    plot(freqs4, spectra4(1,:), 'red-*'); 
    hold on; 
    plot(freqs5, spectra5(1,:), 'green-o'); 
    % Adding labels and title
    xlabel('Frequency (Hz)');
    xlim([0 40]);
    ylabel('Power (dB)');
    title(['Individual Power Spectra during ', condition]);
    % Add a legend to identify which line corresponds to which EEG data
    % legend({'Raw data', 'Detrend', 'Noise', 'ECG', 'EOG'});
    legend({'M1', 'M2', 'A1', 'A2', 'C3'});
    % Enable grid for better visibility
    grid on;
    % Release the hold on the plot
    hold off;

    ith_channel = [indexM1,indexM2,indexA1,indexA2,indexC3]
    figure;
    hold on;
     for j = 1:5
        % Define subplot position and get handle
        ax = subplot(5, 1, j + (1 - 1) * 5);
        % Plot spectrogram for current channel and step
        spectrogram(EEG_raw.data(ith_channel(j),:), hamming(window_length), noverlap, nfft, srate, 'yaxis');
        ylim(fpass);
        caxis([-40, 40]);
        colormap('jet');
        colorbar;
        % Title for each channel
        title(EEG_raw.chanlocs(ith_channel(j)).labels, 'Interpreter', 'none');
     end
     % Set a title for the whole figure
     sgtitle([participant, ' ', condition, ' spectrogram'], 'FontSize', 14);
     hold off;