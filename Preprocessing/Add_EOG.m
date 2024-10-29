  %%% Add the path to your EEGLAB folder
    addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
    eeglab nogui; % Initialize EEGLAB (eeglab nogui; if you want no GUI)

% Define the base folder where the search should start
base_folder = 'D:\Masterarbeit Jannick\Data\GHB_TRA\GHB_TRA_EEG';

% Get a list of all .mat and .edf files in the base folder and its subfolders
all_files = dir(fullfile(base_folder, '**', '*.*'));
all_files = all_files(~[all_files.isdir]); % Filter out directories

% Initialize a cell array to hold grouped files
grouped_files = {};

% Get the folder paths for all files
folder_paths = {all_files.folder};

% Find unique folders where .mat and .edf files are stored
unique_folders = unique(folder_paths);

% Loop through each unique folder
for i = 1:length(unique_folders)
    folder = unique_folders{i};
    
    % Get all .mat and .edf files in this folder
    mat_files_in_folder = dir(fullfile(folder, '*.mat'));
    edf_files_in_folder = dir(fullfile(folder, '*.edf'));
    
    % Find the _preprocessed.mat files
    preprocessed_files = dir(fullfile(folder, '*_preprocessed_EOG.mat'));
    
    % Loop through each _preprocessed.mat file and group with other .mat and .edf files in the same folder
    for j = 1:length(preprocessed_files)
        % Create a new group for the preprocessed file
        group = {};
        group{end+1} = fullfile(preprocessed_files(j).folder, preprocessed_files(j).name); % Add the _preprocessed.mat file
        
        % Add all other .mat files in the folder (if any), excluding the _preprocessed file itself
        other_mat_files = setdiff(fullfile({mat_files_in_folder.folder}, {mat_files_in_folder.name}), group{1});
        if ~isempty(other_mat_files)
            group = [group, other_mat_files]; % Append other .mat files to the group if they exist
        end
        
        % Add any .edf files in the folder (if any)
        if ~isempty(edf_files_in_folder)
            edf_files = fullfile({edf_files_in_folder.folder}, {edf_files_in_folder.name});
            group = [group, edf_files]; % Append .edf files to the group if they exist
        end

        % Add the group to the results
        grouped_files{end+1, 1} = group;
    end
end





% Now process each group by loading the .mat and .edf files
for k = 1:length(grouped_files)
    group = grouped_files{k};
    disp(['Processing Group ' num2str(k)]);

    % Initialize variables
    EEG_clean = [];
    edf_file_loaded = false;  % Flag to check if an EDF file is loaded

    % Loop through the group files
    for file_idx = 1:length(group)
        [~, file_name, ext] = fileparts(group{file_idx}); % Get the file extension

        % Load the single .mat file
        if strcmp(ext, '.mat')
            disp(['Loading MAT file: ' group{file_idx}]);
            load(group{file_idx});  % Load the preprocessed file into EEG_clean
            EEG_clean_file = group{file_idx};  % Store the file name

        % Load .edf files using pop_biosig
        elseif strcmp(ext, '.edf')
            disp(['Loading EDF file: ' group{file_idx}]);
            EEG = pop_biosig(group{file_idx});  % Load the EDF file into EEG_clean
            edf_file_loaded = true;  % Set flag to indicate EDF file is loaded
        end
    end

    % Find the indices for E1 and E2 channels using cellfun
    loc_channel_idx = find(cellfun(@(x) contains(x, 'E1', 'IgnoreCase', true), {EEG.chanlocs.labels}));
    roc_channel_idx = find(cellfun(@(x) contains(x, 'E2', 'IgnoreCase', true), {EEG.chanlocs.labels}));


    % Initialize EEG_clean.data_EOG as an empty array if it isn't already
    if ~isfield(EEG_clean, 'data_EOG')
        EEG_clean.data_EOG = [];
    end

    % Extract LOC and ROC data if available
    if ~isempty(loc_channel_idx) && ~isempty(roc_channel_idx)
        loc_data = EEG.data(loc_channel_idx,:);  % Convert LOC cell to numeric array
        roc_data = EEG.data(roc_channel_idx,:);  % Convert ROC cell to numeric array

        % Combine the LOC and ROC data as two rows (2xX)
        EEG_clean.data_EOG = [loc_data; roc_data];  % LOC as row 1, ROC as row 2
    else
        warning('LOC or ROC channels not found.');
    end

    % Resample the EOG data
    original_srate = EEG.srate;  % Original sampling rate
    target_srate = EEG_clean.srate;  % Desired resampling rate
    EEG_clean.data_EOG = resample(EEG_clean.data_EOG', target_srate, original_srate)';  % Resample

    % Extract the channel labels from the EEG structure
    labels = {EEG.chanlocs.labels};  % Convert to a cell array

    % Find the indices for E1 and E2 channels using cellfun
    loc_channel_idx = find(cellfun(@(x) contains(x, 'E1', 'IgnoreCase', true), labels));
    roc_channel_idx = find(cellfun(@(x) contains(x, 'E2', 'IgnoreCase', true), labels));

    % Assign the labels for EOG channels to EEG_clean
    EEG_clean.channelname_EOG = labels([loc_channel_idx, roc_channel_idx]);

    % Save the updated file
    if ~isempty(EEG_clean_file)
        % Create the new filename by appending _EOG
        [folder, base_name, ~] = fileparts(EEG_clean_file);
        new_file_name = fullfile(folder, [base_name '_EOG.mat']);  % Append '_EOG' to the base name

        % Save the updated EEG_clean structure
        save(new_file_name, 'EEG_clean');
        disp(['Saved file as: ' new_file_name]);
    else
        warning('No _preprocessed.mat file found in the group.');
    end
end













% Now process each group by loading the multiple .mat files
for k = 2:length(grouped_files)
    group = grouped_files{k};
    disp(['Processing Group ' num2str(k)]);

    % Initialize variables
    EEG_clean = [];
    eeg_study = [];

    % Define the base names for loading non-preprocessed files
    non_preprocessed_variable_names = {'eeg_study1', 'eeg_study2', 'eeg_study3', 'eeg_study4'};  % Extend as needed

    % Loop through the group files
    non_preprocessed_count = 1;  % Counter for non-preprocessed files
    for file_idx = 1:length(group)
        [~, file_name, ext] = fileparts(group{file_idx}); % Get the file extension

        % Load .mat files
        if strcmp(ext, '.mat')
            disp(['Loading MAT file: ' group{file_idx}]);

            % Check if the file is the _preprocessed.mat file
            if contains(file_name, '_preprocessed')
                load(group{file_idx});  % Load the preprocessed file into EEG_clean
                EEG_clean_file = group{file_idx};  % Store the file name
            else
                % Dynamically assign the variable name for non-preprocessed files
                if non_preprocessed_count <= length(non_preprocessed_variable_names)
                    eval([non_preprocessed_variable_names{non_preprocessed_count} ' = load(group{file_idx});']);
                    non_preprocessed_count = non_preprocessed_count + 1;
                else
                    warning('More non-preprocessed files than expected.');
                end
            end
        end
    end

    % Find the indices for LOC and ROC channels
    loc_channel_idx = find(contains(eeg_study4.eeg_study.channelnames, 'LOC', 'IgnoreCase', true));
    roc_channel_idx = find(contains(eeg_study4.eeg_study.channelnames, 'ROC', 'IgnoreCase', true));

   
% Initialize EEG_clean.data_EOG as an empty array if it isn't already
if ~isfield(EEG_clean, 'data_EOG')
    EEG_clean.data_EOG = [];
end

% Extract LOC and ROC data if available
if ~isempty(loc_channel_idx) && ~isempty(roc_channel_idx)
    % Convert cell contents to numeric arrays and concatenate correctly
    loc_data = cell2mat(eeg_study3.T(loc_channel_idx));  % Convert LOC cell to numeric array
    roc_data = cell2mat(eeg_study3.T(roc_channel_idx));  % Convert ROC cell to numeric array

    % Ensure the data is a row vector (1xX) for each channel
    loc_data = loc_data(:)';  % Force LOC to be 1xX
    roc_data = roc_data(:)';  % Force ROC to be 1xX

    % Combine the LOC and ROC data as two rows (2xX)
    EEG_clean.data_EOG = [loc_data; roc_data];  % LOC as row 1, ROC as row 2
else
    warning('LOC or ROC channels not found.');
end


    
    EEG_clean.data_EOG = double(EEG_clean.data_EOG);
    % downsample
    % Define original and target sampling rates
    original_srate = eeg_study4.eeg_study.header.channel_frequencies(1);  % Assuming original sampling rate is stored here
    target_srate = EEG_clean.srate;  % Desired resampling rate

    % Resample the EOG data using MATLAB's resample function
    EEG_clean.data_EOG = resample(EEG_clean.data_EOG', target_srate, original_srate)';

    % Also append the channel names to EEG_clean.channelname_EOG
    EEG_clean.channelname_EOG = eeg_study4.eeg_study.channelnames([loc_channel_idx, roc_channel_idx]);

    % Find the _preprocessed file and modify its name
    preprocessed_file_name = EEG_clean_file; % We already stored it earlier

    if ~isempty(preprocessed_file_name)
        % Create the new filename by appending _EOG
        [folder, base_name, ~] = fileparts(preprocessed_file_name);
        new_file_name = fullfile(folder, [base_name '_EOG.mat']);  % Append '_EOG' to the base name

        % Save the updated EEG_clean structure to the new file
        save(new_file_name, 'EEG_clean');
        disp(['Saved file as: ' new_file_name]);
    else
        warning('No _preprocessed.mat file found in the group.');
    end
end


load('D:\Masterarbeit Jannick\Data\SWC\SWC_EEG\VP12\1a_import\SWC_swce212_EEG_preprocessed_EOG.mat')