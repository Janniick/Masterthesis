addpath('D:\Masterarbeit Jannick\scripts\2_Preprocessing\eeglab2024.0');
eeglab nogui;
EEG_mat = pop_biosig('D:\Masterarbeit Jannick\Data\DEX-FX_v2\DEX-FX_v2_EEG\VP102\DEX-FX_102.1.edf');
labels = {EEG_mat.chanlocs.labels};

numChannels = length(EEG_mat.chanlocs);  % Total number of channels
pointsToPlot = 10000;  % Number of points to plot per channel
channelsPerFigure = 12;  % Max channels per figure

for i = 1:numChannels
    if mod(i-1, channelsPerFigure) == 0  % Check if new figure is needed
        figure;  % Create new figure
    end
    subplot(3, 4, mod(i-1, channelsPerFigure) + 1);  % Create subplot (3 rows x 4 cols)
    plot(EEG_mat.data(i, 1:pointsToPlot));  % Plot data for current channel
    title(EEG_mat.chanlocs(i).labels);  % Set the title as the channel label
end

% Plotting the channels
% Finding channel indices
M1_index = find(strcmp(labels, 'M1'));  % Find index for M1
M2_index = find(strcmp(labels, 'M2'));  % Find index for A1
a=EEG_mat.srate*3600;
b=EEG_mat.srate*3600*1.01;
figure;  % Open a new figure window
plot(EEG_mat.data(M1_index, a:b), 'b');
hold on;  % Hold the current plot
plot(EEG_mat.data(M2_index, a:b), 'r');
hold off;  % Release the plot hold
title('Comparison of M1 and M2 Channels');
legend('M1', 'M2');

M2_index = find(strcmp(labels, 'M2'));  % Find index for M2
A2_index = find(strcmp(labels, 'A2'));  % Find index for A2

figure;  % Open a new figure window
plot(EEG_mat.data(M2_index, 1:pointsToPlot), 'b');
hold on;  % Hold the current plot
plot(EEG_mat.data(A2_index, 1:pointsToPlot), 'r');
hold off;  % Release the plot hold
title('Comparison of M2 and A2 Channels');
legend('M2', 'A2');
