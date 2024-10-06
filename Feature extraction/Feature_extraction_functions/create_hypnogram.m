function hypnogram = create_hypnogram(sleep_stages_whole)
    num_epochs = length(sleep_stages_whole);
    hypnogram = zeros(1, num_epochs);
    for epoch = 1:num_epochs
        stage = sleep_stages_whole(epoch);
        if stage == 3 || stage == 4
            hypnogram(epoch) = 3;
        elseif stage >= 6
            hypnogram(epoch) = 6;
        else
            hypnogram(epoch) = stage;
        end
    end
end