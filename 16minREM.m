% other functions (maybe here not needed, but for feature extraction at least they are needed
addpath(genpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2023_mesmart\scripts'));

% preprocessing Funktionen
addpath('O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\matlab_processing_benji\');

% eeglab path
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui; % Initialize EEGLAB (eeglab nogui; if you want no GUI)

% plot Funktion
addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing');

% edf file directory
folderpath = 'D:\Masterarbeit Jannick\Data\DEX-FX_v1\DEX-FX_v1_EEG';
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
    thename = split(edf_name, '.');
    thename = thename{1};
    % Where to save the generated files
    [pathtosave, ~, ~] = fileparts(edf_path);
    % participant number
    % Get the participant number
    [parentPath, ~, ~] = fileparts(edf_path);
    [~, participant, ~] = fileparts(parentPath);

    % Use eeglab biosig to load edf eeg file
    fprintf('---load data\n');
    EEG = pop_biosig(edf_path);
    EEG.name=edf_name;

    % downsample
    fprintf('---downsample\n');
    EEG = pop_resample(EEG, 256);

    % select channels
    channelLabels = {EEG.chanlocs.labels};
    ECG = pop_select(EEG, 'channel', {'ECG'});
    EEG = pop_select(EEG, 'channel', {'F3', 'F4'});
    %EEG = pop_select(EEG, 'channel', {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'});
    newchannelLabels = {EEG.chanlocs.labels};
    
    % Define start and end times in minutes for the REM window
    startTime_REM = 2 * 60 + 45;  % Start at 2 hours and 45 minutes
    endTime_REM = startTime_REM + 16;  % End time 16 minutes later
    % Convert times to seconds
    startSecond_REM = startTime_REM * 60;
    endSecond_REM = endTime_REM * 60;
    % Calculate sample indices
    startIndex_REM = round(startSecond_REM * EEG.srate) + 1;
    endIndex_REM = round(endSecond_REM * EEG.srate);
    % Extract EEG data for the REM window
    selectedEEGData_REM = EEG.data(:, startIndex_REM:endIndex_REM);
    % Create a new EEG structure for the REM window
    EEG_REM = EEG;  % Copy the original structure
    EEG_REM.data = selectedEEGData_REM;  % Replace data with the REM window
    EEG_REM.pnts = size(selectedEEGData_REM, 2);  % Update the number of points
    EEG_REM.xmax = EEG_REM.xmin + (EEG_REM.pnts - 1) / EEG.srate;  % Update xmax accordingly
    % Extract ECG data for the REM window
    selectedECGData = ECG.data(:, startIndex_REM:endIndex_REM);
    % Create a new ECG structure for the selected window
    ECG_REM = ECG;  % Copy the original structure
    ECG_REM.data = selectedECGData;  % Replace data with the selected window
    ECG_REM.pnts = size(selectedECGData, 2);  % Update the number of points
    ECG_REM.xmax = ECG_REM.xmin + (ECG_REM.pnts - 1) / ECG.srate;  % Update xmax accordingly

    % define parameters for plotting
    srate = EEG.srate;
    window_length = 4*srate;
    noverlap = window_length/2;
    nfft = []; % Use default NFFT
    fpass = [0.5 40];
    % how to divide seconds to get values in the plotting frame
    if size(EEG_REM.data(1,:),2) < srate*60 % less than a minute (in seconds)
        divideby = 1;
    elseif size(EEG_REM.data(1,:),2) < srate*60*60 % less than an hour (in minutes)
        divideby = 60;
    else % less than an hour (in hours)
        divideby = 60*60;
    end
    % Number of channels and preprocessing steps
    num_channels = length(EEG_REM.chanlocs);
    num_preprocessing_steps = 6; % Adjust this as needed
    num_preprocessing_steps_array = 1:num_preprocessing_steps;
    % Dynamic positioning of titles in plot
    row_positions = linspace(0.9, 0.1, num_preprocessing_steps+1);

    % start of plotting
    figure('visible', 'on');
    % Hold the current plot to retain all subplots
    hold on;
    
    preprocessing_step = 1;
    % plot raw data in first row
    for ith_channel = 1:num_channels
    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            title(EEG.chanlocs(ith_channel).labels, 'Interpreter', 'none');
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'Raw data of REM sleep', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end

    %%%%% bandpass and detrend, to remove sweat and other artifacts
    fprintf('---bandpass detrend\n');
    preprocessing_step = 2;
    EEG_REM_detrend = EEG_REM;
    if size(EEG_REM.data,1) > 0
        for chcha = 1:size(EEG_REM.data,1)
            EEG_REM_detrend.data(chcha,:) = detrend(EEG_REM.data(chcha,:));
            EEG_REM_detrend.data(chcha,:) = bandpass(EEG_REM.data(chcha,:), [0.5, 40], EEG.srate);
        end
    end
    
    % plot detrended data in second row
    for ith_channel = 1:num_channels
    
    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM_detrend.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'Bandpass & Detrend', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end      

    %%%%% detect bad signal to remove from signal
    preprocessing_step = 3;
    [remove_channel_which, cut_chunks, rem_ind_all] = eeg_acrossfreq_artifact(EEG_REM_detrend);
    
    % plotting after removal of bad signal
    for ith_channel = 1:length(EEG.chanlocs)
    % cut chunks
            scutc = size(cut_chunks{ith_channel},1);
            ccunk = cut_chunks{ith_channel};

    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM_detrend.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'Removal of bad signal', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

    %generate all the black boxes to overlay the bad signal
    for ccc = 1:length(ccunk)
        x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);
    end
    end

   

    %%%%% remove electrical noise from sienna, and make sure to exclude the movement artifacts for calculation (output EEG still has all data, including the movments)
    preprocessing_step = 4;
    EEG_REM_noise = EEG_REM_detrend;
    EEG_REM_noise = eeg_electricalnoise_reducer(EEG_REM_detrend, rem_ind_all);
    % plotting after removal of electrical noise
    for ith_channel = 1:length(EEG.chanlocs)
    % cut chunks
            scutc = size(cut_chunks{ith_channel},1);
            ccunk = cut_chunks{ith_channel};

    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM_noise.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'Removal of electrical noise', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

    %generate all the black boxes to overlay the bad signal
    for ccc = 1:length(ccunk)
        x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);
    end
    end

    

    %%%%% remove heart beat artifacts here again I require the movement artifacts to be excluded (eeg_jelena_ecg_removal will probably not yet run correctly, I will need to adjust it)
    preprocessing_step = 5;
    EEG_REM_ECG = EEG_REM_noise;
    EEG_REM_ECG = eeg_jelena_ecg_removal(EEG_REM_noise, ECG_REM, cut_chunks);
    % plotting after removal of ECG
    for ith_channel = 1:length(EEG.chanlocs)
    % cut chunks
            scutc = size(cut_chunks{ith_channel},1);
            ccunk = cut_chunks{ith_channel};

    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM_ECG.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'ECG removal', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

    %generate all the black boxes to overlay the bad signal
    for ccc = 1:length(ccunk)
        x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);
    end
    end


    %%%%% do some EOG removal here
    preprocessing_step = 6;
    EEG_REM_EOG = EEG_REM_ECG;
    % EEG_REM_EOG = pop_dusk2dawn_clean(EEG_REM_ECG);
    [~,~,EEG_REM_EOG,removed_channels] = clean_artifacts(EEG_REM_EOG, ...
        'Highpass','off', ...
        'ChannelCriterion', 'off', ... % 0 or 0.3 correlation? high density EEG is different
        'FlatlineCriterion', 'off', ...
        'LineNoiseCriterion', 'off', ...
        'BurstCriterion', 30, ... %burst criterion!!
        'WindowCriterion', 0.25);
    
    % plotting after removal of EOG
    for ith_channel = 1:length(EEG.chanlocs)
    % cut chunks
            scutc = size(cut_chunks{ith_channel},1);
            ccunk = cut_chunks{ith_channel};

    subplot(num_preprocessing_steps, num_channels, ith_channel + (num_preprocessing_steps_array(preprocessing_step) - 1)* num_channels);
            spectrogram(EEG_REM_EOG.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            annotation('textbox', [0, row_positions(num_preprocessing_steps_array(preprocessing_step)), 1, 0.05], 'String', 'EOG removal', ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

    %generate all the black boxes to overlay the bad signal
    for ccc = 1:length(ccunk)
        x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);
    end
    end
     % Release hold after plotting all subplots
    hold off;
    % Saving the plot as png in pathtosave
    filename = ['Spectrogram_REM_', participant, '.png'];
    print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');

    % Calculate spectrum for each step
    % For EEG_REM (Raw data)
    [spectra1, freqs1] = spectopo(EEG_REM.data, 0, EEG_REM.srate, 'plot', 'off');
    % For EEG_REM_detrend
    [spectra2, freqs2] = spectopo(EEG_REM_detrend.data, 0, EEG_REM_detrend.srate, 'plot', 'off');
    % For EEG_REM_noise
    [spectra3, freqs3] = spectopo(EEG_REM_noise.data, 0, EEG_REM_noise.srate, 'plot', 'off');
    % For EEG_REM_ECG
    [spectra4, freqs4] = spectopo(EEG_REM_ECG.data, 0, EEG_REM_ECG.srate, 'plot', 'off');
    % For EEG_REM_EOG
    [spectra5, freqs5] = spectopo(EEG_REM_EOG.data, 0, EEG_REM_EOG.srate, 'plot', 'off');
   
    % Create a figure for all spectrograms
    figure;
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
    title('Comparison of EEG Power Spectra in REM');
    % Add a legend to identify which line corresponds to which EEG data
    legend({'Raw data', 'Detrend', 'Noise', 'ECG', 'EOG'});
    % Enable grid for better visibility
    grid on;
    % Release the hold on the plot
    hold off;

    % Saving the plot as png in pathtosave
    filename = ['Powerspectrum_individual_REM_', participant, '.png'];
    print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
    
    % Create a figure for all spectrogram-differences
    figure;
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
    title('Comparison of EEG Power Spectra in REM');
    % Add a legend to identify which line corresponds to which EEG data
    legend({'Raw data - Detrend', 'Detrend - Noise', 'Noise - ECG', 'ECG - EOG'});
    % Enable grid for better visibility
    grid on;
    % Release the hold on the plot
    hold off;

    % Saving the plot as png in pathtosave
    filename = ['Powerspectrum_differences_REM_', participant, '.png'];
    print(gcf, fullfile(pathtosave, filename ), '-dpng', '-r300');
    

    save('EEG_preprocessed.mat', 'EEG'); % later in another feature extraction script you could load('EEG_preprocessed.mat'), of course by including the edf_name in the name

end




    % Define start and end times in minutes for the Non-REM window
    startTime_nonREM = 1 * 60 + 30;  % Start at 1:30 hours
    endTime_nonREM = startTime_nonREM + 16;  % End time 16 minutes later
    % Convert times to seconds
    startSecond_nonREM = startTime_nonREM * 60;
    endSecond_nonREM = endTime_nonREM * 60;
    % Calculate sample indices
    startIndex_nonREM = round(startSecond_nonREM * EEG.srate) + 1;
    endIndex_nonREM = round(endSecond_nonREM * EEG.srate);
    % Extract EEG data for the Non-REM window
    selectedEEGData_nonREM = EEG_REM.data(:, startIndex_nonREM:endIndex_nonREM);
    % Create a new EEG structure for the Non-REM window
    EEG_nonREM = EEG;  % Copy the original structure
    EEG_nonREM.data = selectedEEGData_nonREM;  % Replace data with the Non-REM window
    EEG_nonREM.pnts = size(selectedEEGData_nonREM, 2);  % Update the number of points
    EEG_nonREM.xmax = EEG_nonREM.xmin + (EEG_nonREM.pnts - 1) / EEG.srate;  % Update xmax accordingly