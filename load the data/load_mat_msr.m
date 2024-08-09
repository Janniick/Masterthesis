EEG_mat = eeg_emptyset; % creates empty EEGLAB structure
    ldata = load("D:\Masterarbeit Jannick\Data\PSS\PSS_EEG\VP002\1a_import\pss002_bn_info.mat"); % reads matlab file (save or not??)

    EEG_mat.srate = ldata.eeg_study.header.channel_frequencies(1);
    EEG_mat.data = ldata.eeg_study.data;


    % if data uneven, change to even
    if mod(size(EEG_mat.data,2),2)
        EEG_mat.data = EEG_mat.data(:,1:(end-1));
    end

    EEG_mat.nbchan = length(ldata.eeg_study.channelnames);
    EEG_mat.pnts = size(EEG_mat.data,2);
    EEG_mat.times = single(0:1/EEG_mat.srate:(size(EEG_mat.data, 2) - 1)/ EEG_mat.srate);
    EEG_mat.chanlocs = struct('labels', cellstr(ldata.eeg_study.channelnames));
    EEG_mat = eeg_checkset(EEG_mat);

    ldata.eeg_study.channelnames = ldata.eeg_study.channelnames.';
    