% author : Benjamin Stucky
% year   : 2024
function new_cuts_idx = eeg_cut_and_glue(EEG, cut_chunks, chan, look_around_s)
    arguments
        EEG % and EEGLAB structure
        cut_chunks % a struct with a matrix per channel, with first column from and second column to cuts
        chan (1,1) {mustBeInteger,mustBeReal, mustBeGreaterThanOrEqual(chan, 1)}
		look_around_s  (1,1) {mustBeNumeric,mustBeReal} = 0.8
    end
    % print what one is doing
    %fprintf('---Benjis cut and glue algorithm: based on most similar zero crossings\n');

    new_cuts_idx = zeros(size(cut_chunks{chan},1),2);
    if ~isempty(cut_chunks{chan})
        % loop over all omits
        for gl = 1:size(cut_chunks{chan},1)
            original_cut = cut_chunks{chan}(gl,:);

            % find zero crossings left and right
            left_look = floor(max(original_cut(1)-look_around_s,1)*EEG.srate):floor(original_cut(1)*EEG.srate);
            right_look = floor(max(original_cut(2),1)*EEG.srate):floor(min(original_cut(2)+look_around_s,length(EEG.data(chan,:))/EEG.srate)*EEG.srate);

            zero_crossings_left = find(diff(sign(EEG.data(chan,left_look))));
            zero_crossings_right = find(diff(sign(EEG.data(chan,right_look))));

            % if at the end or no zero crossings, keep cuts
            if isempty(zero_crossings_left)|isempty(zero_crossings_right)|(min(left_look) == 1)|(max(right_look) == length(EEG.data(chan,:)))
                % keep original cuts
                new_cuts_idx(gl,:) = [original_cut(1)*EEG.srate, original_cut(2)*EEG.srate];
            else
                % zero crossing distance
                importance_scale = 0; % 0 means distance will not play a role in minimum, 1 means 1-1 linear,...
                zcl_dist = 1+(1-(left_look-min(left_look))/EEG.srate)*importance_scale;
                zcr_dist = 1+((right_look-min(right_look))/EEG.srate)*importance_scale;


                % derivatives
                grad_1_left = gradient(EEG.data(chan,left_look));
                grad_2_left = gradient(gradient(EEG.data(chan,left_look)));
                grad_3_left = gradient(gradient(gradient(EEG.data(chan,left_look))));
                grad_4_left = gradient(gradient(gradient(gradient(EEG.data(chan,left_look)))));

                grad_1_right = gradient(EEG.data(chan,right_look));
                grad_2_right = gradient(gradient(EEG.data(chan,right_look)));
                grad_3_right = gradient(gradient(gradient(EEG.data(chan,right_look))));
                grad_4_right = gradient(gradient(gradient(gradient(EEG.data(chan,right_look)))));


                % % first right, which min on left hand side?
                % [gmin, imin] = min(abs(grad_1_left(zero_crossings_left)-grad_1_right(zero_crossings_right(1))) + ...
                %     abs(grad_2_left(zero_crossings_left)-grad_2_right(zero_crossings_right(1))) + ...
                %     abs(grad_3_left(zero_crossings_left)-grad_3_right(zero_crossings_right(1))) + ...
                %     abs(grad_4_left(zero_crossings_left)-grad_4_right(zero_crossings_right(1))));
                %
                % % new cuts
                % new_cuts_idx(gl,:) = [left_look(zero_crossings_left(imin)) + 1, right_look(zero_crossings_right(1))]; % + 1 because do not keep two zero crossing points
                %

                % Create 2D arrays for all combinations of zero_crossings_left and zero_crossings_right
                [ZCL, ZCR] = meshgrid(zero_crossings_left, zero_crossings_right);

                % Compute absolute differences for all combinations
                abs_diffs = abs(grad_1_left(ZCL).*zcl_dist(ZCL) - grad_1_right(ZCR).*zcr_dist(ZCR)).^2 + ...
                    abs(grad_2_left(ZCL).*zcl_dist(ZCL) - grad_2_right(ZCR).*zcr_dist(ZCR)).^2 + ...
                    abs(grad_3_left(ZCL).*zcl_dist(ZCL) - grad_3_right(ZCR).*zcr_dist(ZCR)).^2 + ...
                    abs(grad_4_left(ZCL).*zcl_dist(ZCL) - grad_4_right(ZCR).*zcr_dist(ZCR)).^2;

                % Find the minimum value and its indices
                [gmin, imin] = min(abs_diffs(:));

                % Convert linear index to subscript indices
                [row, col] = ind2sub(size(abs_diffs), imin);

                if size(ZCL,2) == 1 % if ax1 matrix, then grad_1_left etc collapses to 1xa vector not ax1 matrix
                    col_o = col;
                    col = row;
                    row = col_o;
                end

                % new cuts
                new_cuts_idx(gl,:) = [left_look(ZCL(row,col)) + 1, right_look(ZCR(row,col))]; % + 1 because do not keep two zero crossing points


                % figure;
                % subplot(2,1,1) % This means 2 rows, 1 column, first plot
                % plot_o_cut = original_cut*EEG.srate - min(left_look);
                % plot(EEG.data(chan,min(left_look):max(right_look)));
                % line([plot_o_cut(1), plot_o_cut(1)], ylim, 'Color', [1 0.5 0]);
                % line([plot_o_cut(2), plot_o_cut(2)], ylim, 'Color', [1 0.5 0]);
                % line([left_look(ZCL(row,col)), left_look(ZCL(row,col))] + 1 - min(left_look), ylim, 'Color', 'r');
                % line([right_look(ZCR(row,col)), right_look(ZCR(row,col))] - min(left_look), ylim, 'Color', 'r');
                % subplot(2,1,2) % This means 2 rows, 1 column, second plot
                % takoutt = new_cuts_idx(gl,1):new_cuts_idx(gl,2);
                % plot(EEG.data(chan,setdiff(min(left_look):max(right_look), takoutt)));
                % line([new_cuts_idx(gl,1) - min(left_look), new_cuts_idx(gl,1) - min(left_look)], ylim, 'Color', 'r');


            end
        end
    else
        new_cuts_idx = [];
    end
    new_cuts_idx(new_cuts_idx == 0) = 1;	
end