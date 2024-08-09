EEG_mat = eeg_emptyset; % creates empty EEGLAB structure
    ldata = load("D:\Masterarbeit Jannick\Data\SWC\SWC_EEG\VP44\1a_import\swce144.mat"); % reads matlab file (save or not??)

    EEG_mat.srate = ldata.outdata.stages.sample_rate_hz(1);
    EEG_mat.data = ldata.outdata.data;


    % if data uneven, change to even
    if mod(size(EEG_mat.data,2),2)
        EEG_mat.data = EEG_mat.data(:,1:(end-1));
    end

    EEG_mat.nbchan = length(ldata.outdata.channel_names);
    EEG_mat.pnts = size(EEG_mat.data,2);
    EEG_mat.times = single(0:1/EEG_mat.srate:(size(EEG_mat.data, 2) - 1)/ EEG_mat.srate);
    EEG_mat.chanlocs = struct('labels', cellstr(ldata.outdata.channel_names));
    EEG_mat = eeg_checkset(EEG_mat);

    clear ldata;