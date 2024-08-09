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

    % Use eeglab biosig to load edf eeg file
    fprintf('---load data\n');
    EEG = pop_biosig(edf_path);
    EEG.name=edf_name;

    % downsample
    fprintf('---downsample\n');
    EEG = pop_resample(EEG, 256);

    % define parameters for plotting
    srate = EEG.srate;
    window_length = 4*srate;
    noverlap = window_length/2;
    nfft = []; % Use default NFFT
    fpass = [0.5 40];
    % how to divide seconds to get values in the plotting frame
    if size(EEG.data(1,:),2) < srate*60 % less than a minute (in seconds)
        divideby = 1;
    elseif size(EEG.data(1,:),2) < srate*60*60 % less than an hour (in minutes)
        divideby = 60;
    else % less than an hour (in hours)
        divideby = 60*60;
    end

    % select channels
    channelLabels = {EEG.chanlocs.labels};
    ECG = pop_select(EEG, 'channel', {'ECG'});
    EEG = pop_select(EEG, 'channel', {'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'});
    newchannelLabels = {EEG.chanlocs.labels};

    % start of plotting
    figure('visible', 'on');
    % Hold the current plot to retain all subplots
    hold on;

    % plot raw data in first row
    for ith_channel = 1:length(EEG.chanlocs)
    % subplot(nr_of_preprocessing_steps, nr_channels, ith_channel + (ith_preprocessing_step - 1)* nr_channels);
    subplot(6, length(EEG.chanlocs), ith_channel + (1 - 1)* length(EEG.chanlocs));
            spectrogram(EEG.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
            title(EEG.chanlocs(ith_channel).labels, 'Interpreter', 'none');
    end

    %%%%% bandpass and detrend, to remove sweat and other artifacts
    fprintf('---bandpass detrend\n');
    if size(EEG.data,1) > 0
        for chcha = 1:size(EEG.data,1)
            EEG.data(chcha,:) = detrend(EEG.data(chcha,:));
            EEG.data(chcha,:) = bandpass(EEG.data(chcha,:), [0.5, 40], EEG.srate);
        end
    end
    
    % plot detrended data in second row
    for ith_channel = 1:length(EEG.chanlocs)
    
     subplot(6, length(EEG.chanlocs), ith_channel + (2 - 1)* length(EEG.chanlocs));
            spectrogram(EEG.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;
    end      

    % Release hold after plotting all subplots
    hold off;

    %%%%% detect bad signal to remove from signal
    [remove_channel_which, cut_chunks, rem_ind_all] = eeg_acrossfreq_artifact(EEG);

    %%%%% remove electrical noise from sienna, and make sure to exclude the movement artifacts for calculation (output EEG still has all data, including the movments)
    EEG = eeg_electricalnoise_reducer(EEG, rem_ind_all);

    %%%%% remove heart beat artifacts here again I require the movement artifacts to be excluded (eeg_jelena_ecg_removal will probably not yet run correctly, I will need to adjust it)
    EEG = eeg_jelena_ecg_removal(EEG, ECG, cut_chunks);

    %%%%% do some EOG removal here

    ...
    prep_steps_plot(EEG);

    save('EEG_preprocessed.mat', 'EEG'); % later in another feature extraction script you could load('EEG_preprocessed.mat'), of course by including the edf_name in the name

end


   for ith_channel = 1:length(EEG.chanlocs)
    % cut chunks
            scutc = size(cut_chunks{ ith_channel},1);
            ccunk = cut_chunks{ith_channel};

    subplot(nr_of_preprocessing, nr_channels, ith_channel + (ith_preprocessing_step - 1)* nr_channels);
            spectrogram(EEG.data(ith_channel,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40, 40]);
            colormap('jet');
            colorbar;

    x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);

   end


 



%%%%% save plot on disc
fig = gcf; % get current figure
        fig.Position(3) = nrchio*600;
        fig.Position(4) = nmat*160;
        exportgraphics(fig, ...
            cfg.savenamepath, ...
            'Resolution',300)