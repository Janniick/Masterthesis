% Author: Benjamin Stucky
% Year: 2024

% Input: give overall path to HDR (and .r* and .vis files) and path to
%        store the resulting data
% Output: Staging, EEG Data, Channel Names

path_to_data = 'O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2024_jannick_sleepquality\ACJ\';
path_for_saving = 'O:\BenjaminS\Benjamin Stucky Projekte und Skripte\Projekte\2024_jannick_sleepquality\';


%addpath 'o:\matlab\' -begin
%addpath 'O:\Diego\myScripts_matlab\ACO-ACJ\Stages'

% list all HDR files
hdr_files = dir(fullfile(path_to_data, '**/*.HDR'));
hdr_id_low = cellfun(@(x)strsplit(x,'.'),lower({hdr_files.name}),'UniformOutput',false);
hdr_id_low = cellfun(@(x)x(1), hdr_id_low);
% list all vis files
vis_files = dir(fullfile(path_to_data, '**/*.vis'));
% list all r* files
r_files = dir(fullfile(path_to_data, '**/*'));
r_files = r_files(~cellfun(@isempty, regexpi({r_files.name}, '\.r\d+$', 'once')));
r_id_low = cellfun(@(x)strsplit(x,'.'),lower({r_files.name}),'UniformOutput',false);
r_id_low = cellfun(@(x)x(1), r_id_low);


% loop over vis files (scoring info), this is because some adaptation night have no
% scoring etc.
for fi = 1:size(vis_files,1)
    % extract info
    vis_name = vis_files(fi).name;
    id_low = lower(vis_name);
    id_low = strsplit(id_low,'.');
    id_low = id_low{1};
    vis_folder = vis_files(fi).folder;



    vis_path = fullfile(vis_folder, vis_name);
    hdr_path = fullfile(hdr_files(strcmp(hdr_id_low, id_low)).folder,hdr_files(strcmp(hdr_id_low, id_low)).name);
    r_path = fullfile(r_files(strcmp(r_id_low, id_low)).folder,r_files(strcmp(r_id_low, id_low)).name);


    hrd_content = splitlines(fileread(hdr_path));

    epoch_len = split(hrd_content(contains(hrd_content, 'Duration')), 'x ');
    first_number_epoch = regexp(epoch_len{1,2}, '\d*', 'match');

    % duration_number = regexp(epoch_len{1,1}, '\d+', 'match');
    duration_number = split(epoch_len{1}, '(');
    duration_number = str2double(duration_number(end)); % number of 20s epochs
    epoch_len = str2num(first_number_epoch{1});

    % channel info
    cha_from = find(contains(hrd_content, '---'));
    cha_from = cha_from(1)+1;
    cha_to = find(strcmp(hrd_content,'.'));
    cha_to = cha_to(1)-1;
    nr_channels = cha_to-cha_from + 1;
    cha_info= cellfun(@strsplit, hrd_content(cha_from:cha_to), 'UniformOutput', false);
    cha_scaling = cellfun(@(x)str2double(x(2)), cha_info);
    cha_names = cellfun(@(x)x(1), cha_info);
    cha_names = strrep(cha_names, '*', '');



    sampling = hrd_content(~cellfun(@isempty, regexp(hrd_content, 'Sampling.*Hz')));
    nsample = regexp(sampling, '\d+', 'match');
    nsample = str2double(nsample{1});

    % read in staging data
    vis_content = splitlines(fileread(vis_path)); % read in vis file, get all rows
    vis_offset = str2double(vis_content{1}); % the offset should be always at first position
    vis_content = vis_content(2:end); % remove first position

    % split each line at spaces
    vis_content = cellfun(@(line) split(line, ' '), vis_content, 'UniformOutput', false);

    % find stopping signal 'e'
    stop_signal = find(cellfun(@(cell) any(strcmp(cell, 'e')), vis_content));
    vis_content = vis_content(1:stop_signal);

    epoch_20s = cellfun(@(cell) cell{1}, vis_content, 'UniformOutput', false);
    stage_20s = cellfun(@(cell) cell{2}, vis_content, 'UniformOutput', false);

    staging_data = table(cellfun(@str2double, epoch_20s), cell2mat(stage_20s), ...
        repelem(epoch_len,size(epoch_20s,1))', repelem(nsample,size(epoch_20s,1))', ...
        'VariableNames', {'epoch', 'stage', 'duration_s', 'sample_rate_hz'});

    staging_data.stage = replace(string(staging_data.stage), 'r', '5'); % rem to 5
    staging_data.stage = replace(staging_data.stage, 'm', '7'); % movement to 7
    staging_data.stage = replace(string(staging_data.stage), '4', '3'); % NREM 4 to 3
    staging_data.stage = replace(string(staging_data.stage), 'e', 'stop'); % NREM 4 to 3



    % get out eeg data
    sampfreq = staging_data.sample_rate_hz(1);
    epochdur = staging_data.duration_s(1);
    datalength = duration_number*epochdur*sampfreq;
    first_ep = 1;

    offset=vis_offset*duration_number*sampfreq*nr_channels*2;           % how many BYTES are to be jumped
    jumpbytes=(first_ep-1)*duration_number*sampfreq*nr_channels*2;      % bytes to go to first_ep

    fidraw = fopen(r_path, 'r');
    fseek(fidraw, offset+jumpbytes, 'bof');
    [unscaleddata, count]=fread(fidraw,[nr_channels, datalength],'short');	% read raw data

    fclose(fidraw);

    % scale the channels
    % two types of scalings (Peter Achermann said it changed at some point, but this calculates it automatically), converts to uV
    if any(cha_scaling>1)
        cha_scaling=1./cha_scaling/2048*1e6; % altes format
    else
        cha_scaling=1e6/2048.*cha_scaling;   % neues format
    end

    % apply the scaling
    for i = (1:nr_channels)
        unscaleddata(i,:) = unscaleddata(i,:)*cha_scaling(i);
    end
    eeg_data = unscaleddata;
    clear unscaleddata;

    eeg_data = single(eeg_data);

    % save data
    outdata = struct('stages', staging_data, 'data', eeg_data, 'channel_names', cha_names);

    save(fullfile(path_for_saving,[id_low,'.mat']),"outdata", 	'-v7.3')
end

%load(fullfile(path_for_saving,[id_low,'.mat']))




