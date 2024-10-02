function ordered_sequence = symb_king(values)
    % Check if the input sequence has at least 3 elements
    if numel(values) < 3
        error('Input sequence must have at least 3 elements.');
    end

    % Initialize the output sequence
    ordered_sequence = [];
    
    % Iterate through overlapping triplets
    for i = 1:3:numel(values)
        triplet = values(i:min(i+2, numel(values)));

        rank_t = tiedrank(triplet);

        % Determine the ordering type
        if all(rank_t == [3, 2, 1])
            ordered_sequence = [ordered_sequence, 1]; %a
        elseif all(rank_t == [3, 1, 2])
            ordered_sequence = [ordered_sequence, 2]; %b
        elseif all(rank_t == [2, 3, 1])
            ordered_sequence = [ordered_sequence, 3]; %c
        elseif all(rank_t == [2, 1, 3])
            ordered_sequence = [ordered_sequence, 4]; %d
        elseif all(rank_t == [1, 2, 3])
            ordered_sequence = [ordered_sequence, 5]; %e 
        else
            ordered_sequence = [ordered_sequence, 6]; %f
        end
    end
end