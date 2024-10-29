 %%% Add the path to your EEGLAB folder
    addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
    eeglab nogui; % Initialize EEGLAB (eeglab nogui; if you want no GUI)

    %%%%% GHB02
load('D:\Masterarbeit Jannick\Data\GHB02\GHB02_EEG\VP07\1a_import\GHB02_vp07_exp1_EEG_preprocessed_EOG.mat');
%test1 = load('D:\Masterarbeit Jannick\Data\GHB02\GHB02_EEG\VP04\1a_import\vp04_exp1_night_datapart_1.mat');
test2 = load('D:\Masterarbeit Jannick\Data\GHB02\GHB02_EEG\VP04\1a_import\vp04_exp1_night_datapart_2.mat');
info = load('D:\Masterarbeit Jannick\Data\GHB02\GHB02_EEG\VP04\1a_import\vp04_exp1_night_info.mat');

EEG_clean.data_ECG = double(test2.T{27})';
original_srate = info.eeg_study.header.channel_frequencies(1);
target_srate = 128;
EEG_clean.data_ECG = resample(EEG_clean.data_ECG, target_srate, original_srate);

save('D:\Masterarbeit Jannick\Data\GHB02\GHB02_EEG\VP04\1a_import\GHB02_vp04_exp1_EEG_preprocessed_EOG.mat', 'EEG_clean');



%%%%GHB_TRA
filepath = 'D:\Masterarbeit Jannick\Data\GHB_TRA\GHB_TRA_EEG\H19\GHBTRA_H19_EEG_preprocessed_EOG.mat';
load(filepath);
EEG_raw = pop_biosig('D:\Masterarbeit Jannick\Data\GHB_TRA\GHB_TRA_EEG\H19\MGT_H19_N5_sleep_sienna3_749r0s1.edf')
ECG_idx = find(cellfun(@(x) contains(x, 'ECG', 'IgnoreCase', true), {EEG_raw.chanlocs.labels}));
EEG_clean.data_ECG = EEG_raw.data(ECG_idx,:);
original_srate = EEG_raw.srate;
target_srate = 128;
EEG_clean.data_ECG = resample(EEG_clean.data_ECG, target_srate, original_srate);
save(filepath, 'EEG_clean');