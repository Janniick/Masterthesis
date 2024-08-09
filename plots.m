function prep_steps_plot(varargin)
    % inputs need to be any number >=1 of EEGlab structures with same channels and sampling rate to plot, 
    % at the end with a cfg structure

    % Number of matrices
    nmat = nargin - 1;

    % The last input is the config file
    cfg = varargin{end};
    
    % check cfg, or set defaults
    if ~isfield(cfg,'maintitle') % display the plot
    	cfg.maintitle = 'Preprocessing Steps';
    end
    if ~isfield(cfg,'subtitles') % display the plot
    	cfg.subtitles = arrayfun(@(x) sprintf('step_%d', x), 1:nmat, 'UniformOutput', false);
    end
    if ~isfield(cfg,'displayit') % display the plot
    	cfg.displayit = false;
    end
    if ~isfield(cfg,'saveit') % save the plot
    	cfg.saveit = false;
    end
    if ~isfield(cfg,'savenamepath') % saving name and path
    	cfg.savenamepath = 'prep_steps.png';
    end
    if ~isfield(cfg,'omit_lines_idx_EEGs') % which EEGs to show some omits
    	cfg.omit_lines_idx_EEGs = [];
    end
    if ~isfield(cfg,'omit_lines_cutchunks') % the omits in black lines
    	cfg.omit_lines_cutchunks = [];
    end




    %%%%%% start plotting
    % setup
    if cfg.displayit
        figure('visible', 'on');
    else
        figure('visible', 'off');
    end

    srate = varargin{1}.srate;
    window_length = 4*srate;
    noverlap = window_length/2;
    nfft = []; % Use default NFFT
    fpass = [0.5 40];

    sgtitle(cfg.maintitle);
    nrchio = length(varargin{1}.chanlocs);

    % how to divide seconds to get values in the plotting frame
    if size(varargin{1}.data(1,:),2) < srate*60 % less than a minute (in seconds)
        divideby = 1;
    elseif size(varargin{1}.data(1,:),2) < srate*60*60 % less than an hour (in minutes)
        divideby = 60;
    else % less than an hour (in hours)
        divideby = 60*60;
    end

    for io =  1:nrchio % loop over channels
        channelname = varargin{1}.chanlocs(io).labels;

        if ~(isempty(cfg.omit_lines_cutchunks)|isempty(cfg.omit_lines_idx_EEGs))
            % cut chunks
            scutc = size(cfg.omit_lines_cutchunks{io},1);
            ccunk = cfg.omit_lines_cutchunks{io};
        end

        for isu = 1:nmat % loop over EEGs
            subplot(nmat, nrchio, io + (isu - 1)*nrchio);
            spectrogram(varargin{isu}.data(io,:), hamming(window_length), noverlap, nfft, srate, 'yaxis')
            ylim(fpass);
            caxis([-40 40]);
            colormap('jet');
            colorbar;
            if isu == 1
                title([{['\fontsize{14}', channelname], ['\fontsize{9}',cfg.subtitles{isu}]}]);
            else
                title(cfg.subtitles{isu})
            end

            % plot horizontal lines if indicated
            if ismember(isu,cfg.omit_lines_idx_EEGs)
                % ax = gca;
                for  ccc = 1:scutc
                    x = [ccunk(ccc, 1)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 2)/divideby, ccunk(ccc, 1)/divideby];
                    % y = [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)];
                    y = [fpass(1), fpass(1), fpass(2), fpass(2)];
                    patch(x, y, 'black', 'FaceAlpha', 1);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Save the plot as a PNG file
    if cfg.saveit
        fig = gcf;
        fig.Position(3) = nrchio*600;
        fig.Position(4) = nmat*160;
        exportgraphics(fig, ...
            cfg.savenamepath, ...
            'Resolution',300)
    end
end
