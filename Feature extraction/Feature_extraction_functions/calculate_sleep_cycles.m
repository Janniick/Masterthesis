function cycle_table = calculate_sleep_cycles(hypnogram_py)
    % calculate_sleep_cycles - Calculate sleep cycles from a hypnogram using YASA.
    %
    % Syntax: cycle_table = calculate_sleep_cycles(hypnogram_py)
    %
    % Inputs:
    %    hypnogram_py - Python NumPy array of sleep stages per epoch.
    %
    % Outputs:
    %    cycle_table - MATLAB table containing sleep cycle information:
    %                 n_cycle, start_epoch, end_epoch, duration_min.
    %
    % Example:
    %    hypnogram_py = py.numpy.array([0,1,2,2,3,5,2,3,5,2]);
    %    cycle_table = calculate_sleep_cycles(hypnogram_py);
    %
    % Author:
    %    Jannick (Adapted by ChatGPT)
    %
    % Date:
    %    October 2024

    %% Step 1: Ensure YASA is Imported
    if count(py.sys.path, '') == 0
        insert(py.sys.path, int32(0), '');
    end
    try
        py.importlib.import_module('yasa');
    catch ME
        error('Failed to import YASA. Ensure YASA is installed in the Python environment used by MATLAB.\nError Message: %s', ME.message);
    end

    %% Step 2: Find Sleep Periods Using YASA's hypno_find_periods
    try
        sleep_cycles = py.yasa.hypno_find_periods(hypnogram_py, ...
            pyargs('sf_hypno', py.float(1/30), 'equal_length', py.False));
    catch ME
        error('Error in YASA hypno_find_periods: %s', ME.message);
    end

    %% Step 3: Extract Values from the sleep_cycles DataFrame
    try
       values = double(sleep_cycles.values);
    catch ME
        error('Error extracting values from sleep_cycles: %s', ME.message);
    end

    %% Step 4: Collapse Consecutive Stages
    if isempty(values)
        warning('No sleep cycles detected in the hypnogram.');
        cycle_table = table();
        return;
    end

    collapsed_values = [];
    current_stage = values(1, 1);
    current_start = values(1, 2);
    current_length = values(1, 3);

    for idx = 2:size(values, 1)
        if values(idx, 1) == current_stage
            current_length = current_length + values(idx, 3);
        else
            collapsed_values = [collapsed_values; current_stage, current_start, current_length];
            current_stage = values(idx, 1);
            current_start = values(idx, 2);
            current_length = values(idx, 3);
        end
    end
    % Add the last stage
    collapsed_values = [collapsed_values; current_stage, current_start, current_length];

    %% Step 5: Identify Sleep Cycles Based on Stage Transitions
    sleep_cycles_struct = struct('start_epoch', {}, 'end_epoch', {}, 'duration_epochs', {});

    start_cycle = 0;
    end_cycle = 0;

    for idx = 1:size(collapsed_values, 1)
        stage = collapsed_values(idx, 1);
        start_epoch = collapsed_values(idx, 2);
        length_epoch = collapsed_values(idx, 3);

        if stage == 5  % REM stage
            if start_cycle > 0
                if idx < size(collapsed_values, 1) && (collapsed_values(idx + 1, 1) == 2 || collapsed_values(idx + 1, 1) == 3)
                    end_cycle = collapsed_values(idx + 1, 2) - 1;  % End epoch before next NREM
                    duration_in_epochs = end_cycle - start_cycle;
                    sleep_cycles_struct(end+1).start_epoch = start_cycle;
                    sleep_cycles_struct(end).end_epoch = end_cycle;
                    sleep_cycles_struct(end).duration_epochs = duration_in_epochs;
                    start_cycle = 0;  % Reset for next cycle
                end
            end
        elseif stage == 2 || stage == 3  % NREM stages
            if start_cycle == 0
                start_cycle = start_epoch;
            end
        end
    end

    %% Step 6: Convert sleep_cycles_struct to MATLAB Array
    if isempty(sleep_cycles_struct)
        warning('No complete sleep cycles identified.');
        cycle_table = table();
        return;
    end

    sleep_cycles_array = [];
    for idx = 1:length(sleep_cycles_struct)
        sleep_cycles_array = [sleep_cycles_array; sleep_cycles_struct(idx).start_epoch, ...
                              sleep_cycles_struct(idx).end_epoch, ...
                              sleep_cycles_struct(idx).duration_epochs];
    end

    %% Step 7: Convert Duration from Epochs to Minutes
    duration_in_minutes = (sleep_cycles_array(:, 3) * 30) / 60;  % Convert to minutes

    %% Step 8: Create the Cycle Table
    n_cycles = size(sleep_cycles_array, 1);
    cycle_numbers = (1:n_cycles)';
    start_epochs = sleep_cycles_array(:, 1);
    end_epochs = sleep_cycles_array(:, 2);
    duration_minutes = duration_in_minutes;

    cycle_table = table(cycle_numbers, start_epochs, end_epochs, duration_minutes, ...
        'VariableNames', {'n_cycle', 'start_epoch', 'end_epoch', 'duration_min'});

end
