function varargout=natChord_analyze(what, varargin)

addpath('functions/')

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);

project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

switch (what)
    case 'subject_routine'
        % handling input arguments:
        subject_name = 'subj01';
        smoothing_win_length = 25;
        lpf = 0;                            % wether to do lowpass filtering on the EMG data.
        sampling_option = 'whole_sampled';  % sampling option for the natural EMG data
        natural_window_size = 200;          % wn size for natural EMGs sampling in ms.
        wn_spacing = 2;                     % the spacing between windows for the whole_samlped sampling option (i.e. how many windows to skip)
        natural_window_type = 'Rect';       % window shape for the natural EMG sampling. 'Rect' or 'Gaussian'
        vararginoptions(varargin,{'subject_name','smoothing_win_length','lpf','Fpass_lpf','Fstop_lpf', ...
                          'sampling_option','natural_window_size','natural_window_type','wn_spacing'});
        
        % if a cell containing multiple subjects was given:
        if (iscell(subject_name))
            for i = 1:length(subject_name)
                natChord_subj(subject_name{i},'smoothing_win_length',smoothing_win_length,'lpf',lpf, ...
                             'sampling_option',sampling_option,'natural_window_size',natural_window_size, ...
                             'natural_window_type',natural_window_type,'wn_spacing',wn_spacing);
            end
        % if a single subject as a char was given:
        else
            natChord_subj(subject_name,'smoothing_win_length',smoothing_win_length,'lpf',lpf, ...
                         'sampling_option',sampling_option,'natural_window_size',natural_window_size, ...
                         'natural_window_type',natural_window_type,'wn_spacing',wn_spacing);
        end
    
    case 'make_analysis_data'
        % Calculate RT, MT, Mean Deviation for each trial of each subejct
        % and create a struct without the mov signals and save it as a
        % single struct called efc1_all.mat
        
        % getting subject files:
        files = dir(fullfile(project_path, 'analysis', 'natChord_*_raw.tsv'));
        movFiles = dir(fullfile(project_path, 'analysis', 'natChord_*_mov.mat'));
        
        % container to hold all subjects' data:
        ANA = [];
        
        % looping through subjects' data:
        for i = 1:length({files(:).name})
            % load subject data:
            tmp_data = dload(fullfile(files(i).folder, files(i).name));
            tmp_mov = load(fullfile(movFiles(i).folder, movFiles(i).name));
            tmp_mov = tmp_mov.MOV_struct;
            
            mean_dev_tmp = zeros(length(tmp_data.BN),1);
            rt_tmp = zeros(length(tmp_data.BN),1);
            mt_tmp = zeros(size(rt_tmp));
            first_finger_tmp = zeros(size(rt_tmp));
            % loop through trials:
            for j = 1:length(tmp_data.BN)
                % if trial was correct:
                if (tmp_data.trialCorr(j) == 1)
                    % calculate and store mean dev:
                    mean_dev_tmp(j) = calculate_mean_dev(tmp_mov{j}, tmp_data.chordID(j), ...
                                                         tmp_data.baselineTopThresh(j), tmp_data.RT(j), ...
                                                         tmp_data.fGain1(j), tmp_data.fGain2(j), tmp_data.fGain3(j), ...
                                                         tmp_data.fGain4(j), tmp_data.fGain5(j));
                    % calculate and stor rt and mt:
                    [rt_tmp(j),mt_tmp(j),first_finger_tmp(j)] = calculate_rt_mt(tmp_mov{j}, tmp_data.chordID(j), ...
                                                                tmp_data.baselineTopThresh(j), tmp_data.RT(j), ...
                                                                tmp_data.fGain1(j), tmp_data.fGain2(j), tmp_data.fGain3(j), ...
                                                                tmp_data.fGain4(j), tmp_data.fGain5(j));
                
                % if trial was incorrect:
                else
                    % mean dev:
                    mean_dev_tmp(j) = -1;
                    rt_tmp(j) = -1;
                    mt_tmp(j) = -1;
                    first_finger_tmp(j) = -1;
                end
            end
            
            % removing unnecessary fields:
            tmp_data = rmfield(tmp_data,'RT');
            tmp_data = rmfield(tmp_data,'trialPoint');

            % adding the calculated parameters to the subject struct:
            tmp_data.RT = rt_tmp;
            tmp_data.MT = mt_tmp;
            tmp_data.first_finger = first_finger_tmp;
            tmp_data.mean_dev = mean_dev_tmp;
            
            % adding subject data to ANA:
            ANA=addstruct(ANA,tmp_data,'row','force');
        end

        dsave(fullfile(project_path,'analysis','natChord_all.tsv'),ANA);

    case 'make_natural_dist'
        subject_name = 'subj01';
        fs_emg = 2148.1481;         % EMG sampling rate in Hz  
        lpf = 0;                    % flag to do lowpass filtering;
        Fpass_lpf = 20;
        Fstop_lpf = 30;
        natural_window_size = 200;      % window size to sample natural EMG
        sampling_option = 'whole_sampled';      % sampling option to select windows from natural EMGs.
        natural_window_type = 'Rect';   % sampling window type for natural EMGs.
        wn_spacing = 2;                 % sampling spacing for the 'whole_sampled' option.
        vararginoptions(varargin,{'subject_name','lpf','Fpass_lpf','Fstop_lpf', ...
                                  'sampling_option','natural_window_size','natural_window_type','wn_spacing'});
        
        % EMG filters:
        hd = emg_filter_designer(fs_emg, 'filter_type', 'bandpass');   % creates my default bpf (20,500)Hz
        if (lpf == 1)
            hd_lpf = emg_filter_designer(fs_emg, 'filter_type', 'lowpass','Fpass_lpf',Fpass_lpf,'Fstop_lpf',Fstop_lpf);    % creates lpf
        else
            hd_lpf = [];
        end
        
                
        % if a cell containing multiple subjects was given:
        if (iscell(subject_name))
            for i = 1:length(subject_name)        
                % Preprocessing and dealing with the natural EMGs:
                fprintf("Processing natural EMG data...\n\n")
                make_natural_emg(subject_name{i},fs_emg,hd,hd_lpf,natural_window_type,natural_window_size,sampling_option,wn_spacing);
            end
        % if a single subject as a char was given:
        else
            % Preprocessing and dealing with the natural EMGs:
            fprintf("Processing natural EMG data...\n\n")
            make_natural_emg(subject_name,fs_emg,hd,hd_lpf,natural_window_type,natural_window_size,sampling_option,wn_spacing);
        end

    case 'natural_emg_autocorr'
        subject_name = 'subj01';
        wn_size = 5; % win size for the natural data in ms.
        vararginoptions(varargin,{'subject_name','wn_size'})

        sess = {'sess01','sess02'};
        
        % load whole emg:
        emg_nat = load(fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_whole.mat']));
        emg_nat = emg_nat.emg_natural_dist;
        
        figure;
        for i = 1:length(sess)
            a_avg = 0;
            for ch = 1:size(emg_nat{i},2)
                [a,lags] = autocorr(emg_nat{i}(:,ch),NumLags=1000); 
                a_avg = a_avg + a/size(emg_nat{i},2);

                subplot(2,1,i)
                plot(lags*wn_size,a,'LineWidth',0.1)
                hold on
                % ylim([0,1])
                xlabel('lag (ms)')
                ylabel('corr')
                title(sess{i})
            end
            plot(lags*wn_size,a_avg,'k','LineWidth',2);
            hold on;
            yline(0)
            idx = findNearest(a_avg,0.2);
            plot([idx*wn_size idx*wn_size],[0 a_avg(idx)],':k','LineWidth',2)
            plot([0 idx*wn_size],[a_avg(idx) a_avg(idx)],':k','LineWidth',2)
        end


    case 'avg_chord_patterns'
        % handling input arguments:
        subject_name = 'subj01';    % subject ID
        plot_option = 1;                   % flag to plot the avg chord patterns.
        normalize_channels = 0;     % flag to whether normalize the channels by their norms or not.
        vararginoptions(varargin,{'subject_name','plot_option','normalize_channels'});

        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        data = getrow(data, data.sn == str2double(subject_name(end-1:end)));
        
        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};
        
        % containers:
        chord_emg_mat = cell(length(sess),1);
        
        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_name);

        % looping through sessions:
        if (plot_option)
            figure;
        end
        for i = 1:length(sess)
            % get chords of the sessions:
            sess_rows = data.BN>=sess_blocks{i}(1) & data.BN<=sess_blocks{i}(end) & data.sn==str2double(subject_name(end-1:end));
            chords = unique(data.chordID(sess_rows));
            
            % sorting chords in an arbitrary way:
            chords_sorted = [19999, 91999, 99199, 99919, 99991, 29999, 92999, 99299, 99929, 99992];
            [ind,~] = find(chords == chords_sorted);
            chords = [chords_sorted'; chords(setdiff(1:length(chords),ind))];

            % initialize chord_emg_mat
            chord_emg_mat{i} = zeros(length(chords),10);

            % looping through chords of the sessions:
            for j = 1:length(chords)
                % selecting rows of the chord that were correct trials:
                tmp_row = sess_rows & data.chordID==chords(j) & data.trialCorr==1;

                % getting the avg of the emg_hold_avg
                chord_emg_mat{i}(j,1) = mean(data.emg_hold_avg_e1(tmp_row));
                chord_emg_mat{i}(j,2) = mean(data.emg_hold_avg_e2(tmp_row));
                chord_emg_mat{i}(j,3) = mean(data.emg_hold_avg_e3(tmp_row));
                chord_emg_mat{i}(j,4) = mean(data.emg_hold_avg_e4(tmp_row));
                chord_emg_mat{i}(j,5) = mean(data.emg_hold_avg_e5(tmp_row));
                chord_emg_mat{i}(j,6) = mean(data.emg_hold_avg_f1(tmp_row));
                chord_emg_mat{i}(j,7) = mean(data.emg_hold_avg_f2(tmp_row));
                chord_emg_mat{i}(j,8) = mean(data.emg_hold_avg_f3(tmp_row));
                chord_emg_mat{i}(j,9) = mean(data.emg_hold_avg_f4(tmp_row));
                chord_emg_mat{i}(j,10) = mean(data.emg_hold_avg_f5(tmp_row));

                if (normalize_channels)
                    chord_emg_mat{i}(j,:) = chord_emg_mat{i}(j,:) ./ scales(:,i)';
                end
            end
            
            if (plot_option)
                % plotting the chord EMG patterns:
                subplot(1,2,i)
                pcolor([[chord_emg_mat{i}, zeros(size(chord_emg_mat{i},1),1)] ; zeros(1,size(chord_emg_mat{i},2)+1)])
                colorbar
                
                % plot settings:
                ax = gca;
                
                set(ax,'YTick',(1:size(chord_emg_mat{i},1))+0.5)
                set(ax,'YTickLabel',chords)
                
                set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
                set(ax,'XTickLabel',emg_locs_names)
                
                set(gca,'YDir','reverse')
                title(sess{i})
            end
        end

        fprintf("Corr of pattern matrices across sessions = %.4f\n",corr2(chord_emg_mat{1},chord_emg_mat{2}))
        varargout{1} = chord_emg_mat;
        varargout{2} = chords;

    case 'chord_pattern_corr'
        plot_option = 1;
        vararginoptions(varargin,{'plot_option'})

        % getting subject numbers:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        subjects = unique(data.sn);
        
        % tranforming subject numbers to subject names:
        subjects = strcat('subj',num2str(subjects,'%02.f'));

        % defining sessions:
        sess = {'sess01','sess02'};

        vectorized_patterns = [];
        % looping through subjects:
        for i = 1:size(subjects,1)
            disp(subjects(i,:))
            % getting avg chord patterns of subject i:
            [chord_emg_mat,~] = natChord_analyze('avg_chord_patterns','subject_name',subjects(i,:),'plot_option',0,'normalize_channels',1);

            % vectorizing sujbect emg patterns and concatenating:
            tmp = cellfun(@(x) reshape(x',1,[])', chord_emg_mat, 'UniformOutput', false);

            % concatenating sessions and subjects:
            vectorized_patterns = [vectorized_patterns [tmp{:}]];
        end
        
        % single finger:
        sf_patterns = vectorized_patterns(1:10*size(chord_emg_mat{1},2),:);

        % multi finger:
        mf_patterns = vectorized_patterns(10*size(chord_emg_mat{1},2)+1:end,:);

        % patterns correlations:
        corr_sf = triu(corr(sf_patterns),1);
        corr_mf = triu(corr(mf_patterns),1);
        
        % indices of within subject correlations:
        idx = 1:2:100;
        idx = idx(1:size(subjects,1));
        
        corr_within = [];
        corr_across = [];
        % looping through subjects:
        for i = 1:length(idx)
            % within subject correlations:
            % single finger:
            corr_within(i).sf_within =  corr_sf(idx(i),idx(i)+1);

            % multi finger:
            corr_within(i).mf_within = corr_mf(idx(i),idx(i)+1);
            
            % removing within subject
            corr_sf(idx(i),idx(i)+1) = 0;
            corr_mf(idx(i),idx(i)+1) = 0;
        end

        % across subject correlations mean:
        % single finger:
        corr_across(1).sf_across_avg = mean(corr_sf(corr_sf~=0));
        corr_across(1).sf_across_sem = std(corr_sf(corr_sf~=0))/length(corr_sf(corr_sf~=0));

        % multi finger:
        corr_across(1).mf_across_avg = mean(corr_mf(corr_mf~=0));
        corr_across(1).mf_across_sem = std(corr_mf(corr_mf~=0))/length(corr_mf(corr_mf~=0));

        if plot_option
            X = {'within subjects','across subjects'};
            Y = [mean([corr_within.sf_within]),  corr_across.sf_across_avg ; mean([corr_within.mf_within]), corr_across.mf_across_avg];
            errorplus = [std([corr_within.sf_within])/length([corr_within.sf_within]), corr_across.sf_across_sem ; std([corr_within.mf_within])/length([corr_within.mf_within]), corr_across.mf_across_sem];
            
            figure;
            bar_SEM(Y,errorplus,'type','dashplot','xtick_labels',X,'line_groups',0)
            % set(gca,'XTickLabel',X)
        end

        varargout{1} = vectorized_patterns;
        varargout{2} = corr_within;
        varargout{3} = corr_across;
    
    case 'behavior_reliability'
        plot_option = 1;
        vararginoptions(varargin,{'plot_option'})

        % getting subject numbers:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        subjects = unique(data.sn);

        % defining sessions:
        sess = {'sess01','sess02'};

        cat_MD = [];
        % looping through subjects:
        for sn = 1:size(subjects,1)    
            % getting avg mean deviation of chords:
            chords_mean_dev = [];
            for j = 1:length(sess)
                chords = unique(data.chordID(sess_rows));
            
                % sorting chords in an arbitrary way:
                chords_sorted = [19999, 91999, 99199, 99919, 99991, 29999, 92999, 99299, 99929, 99992];
                [ind,~] = find(chords == chords_sorted);
                chords = [chords_sorted'; chords(setdiff(1:length(chords),ind))];
                for i = 1:length(chords)
                    row = data.sn==subjects(sn) & data.BN>=sess_blocks{j}(1) & data.BN<=sess_blocks{j}(end) & data.trialCorr==1 & data.chordID==chords(i);
                    tmp_mean_dev(i) = mean(data.mean_dev(row));
                    
                end
                chords_mean_dev(:,j) = tmp_mean_dev;
            end
        
            % getting avg chord patterns of subject i:
            [chord_emg_mat,~] = natChord_analyze('avg_chord_patterns','subject_name',subjects(i,:),'plot_option',0,'normalize_channels',1);

            % vectorizing sujbect emg patterns and concatenating:
            tmp = cellfun(@(x) reshape(x',1,[])', chord_emg_mat, 'UniformOutput', false);

            % concatenating sessions and subjects:
            vectorized_patterns = [vectorized_patterns [tmp{:}]];
        end
        
        % single finger:
        sf_patterns = vectorized_patterns(1:10*size(chord_emg_mat{1},2),:);

        % multi finger:
        mf_patterns = vectorized_patterns(10*size(chord_emg_mat{1},2)+1:end,:);

        % patterns correlations:
        corr_sf = triu(corr(sf_patterns),1);
        corr_mf = triu(corr(mf_patterns),1);
        
        % indices of within subject correlations:
        idx = 1:2:100;
        idx = idx(1:size(subjects,1));
        
        corr_within = [];
        corr_across = [];
        % looping through subjects:
        for i = 1:length(idx)
            % within subject correlations:
            % single finger:
            corr_within(i).sf_within =  corr_sf(idx(i),idx(i)+1);

            % multi finger:
            corr_within(i).mf_within = corr_mf(idx(i),idx(i)+1);
            
            % removing within subject
            corr_sf(idx(i),idx(i)+1) = 0;
            corr_mf(idx(i),idx(i)+1) = 0;
        end

        % across subject correlations mean:
        % single finger:
        corr_across(1).sf_across_avg = mean(corr_sf(corr_sf~=0));
        corr_across(1).sf_across_sem = std(corr_sf(corr_sf~=0))/length(corr_sf(corr_sf~=0));

        % multi finger:
        corr_across(1).mf_across_avg = mean(corr_mf(corr_mf~=0));
        corr_across(1).mf_across_sem = std(corr_mf(corr_mf~=0))/length(corr_mf(corr_mf~=0));

        if plot_option
            X = {'within subjects','across subjects'};
            Y = [mean([corr_within.sf_within]),  corr_across.sf_across_avg ; mean([corr_within.mf_within]), corr_across.mf_across_avg];
            errorplus = [std([corr_within.sf_within])/length([corr_within.sf_within]), corr_across.sf_across_sem ; std([corr_within.mf_within])/length([corr_within.mf_within]), corr_across.mf_across_sem];
            
            figure;
            bar_SEM(Y,errorplus,'type','dashplot','xtick_labels',X,'line_groups',0)
            % set(gca,'XTickLabel',X)
        end

        varargout{1} = vectorized_patterns;
        varargout{2} = corr_within;
        varargout{3} = corr_across;

    case 'visualize_natural_emg'
        % handling input arguments:
        sampling_option = 'whole_sampled';
        subject_name = 'subj01';
        normalize_channels = 1;             % flag to whether normalize the channels by their norms or not.
        dimensions = [];                    % dimensions of the natural data to show. by default random dimensions are selected.
        vararginoptions(varargin,{'subject_name','sampling_option','normalize_channels','dimensions'});

        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};

        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_name);
        
        if (normalize_channels)
            % normalizing the natural EMGs:
            for i = 1:length(sess)
                emg_dist{i} = emg_dist{i} ./ scales(:,i)';
            end
        end

        % loading subject data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        
        % calculating avg chord patterns:
        [chord_emg_mat, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name,'plot_option',0,'normalize_channels',normalize_channels);
        
        % getting avg mean deviation of chords:
        chords_mean_dev = [];
        c = [];
        for j = 1:length(sess)
            for i = 1:length(chords)
                row = data.BN>=sess_blocks{j}(1) & data.BN<=sess_blocks{j}(end) & data.trialCorr==1 & data.chordID==chords(i) & data.sn==str2double(subject_name(end-1:end));
                tmp_mean_dev(i) = mean(data.mean_dev(row));
            end
            chords_mean_dev(:,j) = tmp_mean_dev;
        end
        
        % emg locations:
        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];
        
        % select 3 random dimensions:
        dims = randperm(size(chord_emg_mat{1},2));
        dims = dims(1:3);
        
        % if user input dimensions:
        if (~isempty(dimensions))
            if length(dimensions) ~= 3
                warning('Input dimensions length must be 3 -> Changed to random dimensions...')
            else
                dims = dimensions;
            end
        end

        % loop on sessions:
        figure;
        for i = 1:length(sess)
            subplot(1,2,i)
            
            % scatter 3D natural EMG dist:
            scatter3(emg_dist{i}(:,dims(1)), emg_dist{i}(:,dims(2)), emg_dist{i}(:,dims(3)), 10, 'filled', 'MarkerFaceColor', [0.6,0.6,0.6], 'HandleVisibility','off');
            xlabel(emg_locs_names(dims(1)))
            ylabel(emg_locs_names(dims(2)))
            zlabel(emg_locs_names(dims(3)))
            title(sess{i})

            hold all;
            
            % mapping mean devs to colormap:
            c = map2color(chords_mean_dev(:,i), autumn);

            % put avg chord patterns on the plot
            for j = 1:size(chord_emg_mat{i},1)
                % in case of single finger chords:
                if (j <= 10)
                    scatter3(chord_emg_mat{i}(j,dims(1)), chord_emg_mat{i}(j,dims(2)), chord_emg_mat{i}(j,dims(3)), 100, 'k', 'filled', 'HandleVisibility','off')
                % in case of multi finger chords:
                else
                    scatter3(chord_emg_mat{i}(j,dims(1)), chord_emg_mat{i}(j,dims(2)), chord_emg_mat{i}(j,dims(3)), 100, 'filled', 'MarkerFaceColor', c(j,:))
                end
            end
            legend(num2str(chords(11:end)))

        end
        colorbar;

    case 'get_scale_factor_emg'
        subject_name = 'subj01';
        vararginoptions(varargin,{'subject_name'});

        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        data = getrow(data, data.sn == str2double(subject_name(end-1:end)));
        
        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};
        
        % containers:
        chord_emg_mat = cell(length(sess),1);

        % looping through sessions:
        for i = 1:length(sess)
            % get chords of the sessions:
            sess_rows = data.BN>=sess_blocks{i}(1) & data.BN<=sess_blocks{i}(end);
            chords = unique(data.chordID(sess_rows));
            
            % sorting chords in an arbitrary way:
            chords_sorted = [19999, 91999, 99199, 99919, 99991, 29999, 92999, 99299, 99929, 99992];
            [ind,~] = find(chords == chords_sorted);
            chords = [chords_sorted'; chords(setdiff(1:length(chords),ind))];

            % initialize chord_emg_mat
            chord_emg_mat{i} = zeros(length(chords),10);

            % looping through chords of the sessions:
            for j = 1:length(chords)
                % selecting rows of the chord that were correct trials:
                tmp_row = sess_rows & data.chordID==chords(j) & data.trialCorr==1;

                % getting the avg of the emg_hold_avg
                chord_emg_mat{i}(j,1) = mean(data.emg_hold_avg_e1(tmp_row));
                chord_emg_mat{i}(j,2) = mean(data.emg_hold_avg_e2(tmp_row));
                chord_emg_mat{i}(j,3) = mean(data.emg_hold_avg_e3(tmp_row));
                chord_emg_mat{i}(j,4) = mean(data.emg_hold_avg_e4(tmp_row));
                chord_emg_mat{i}(j,5) = mean(data.emg_hold_avg_e5(tmp_row));
                chord_emg_mat{i}(j,6) = mean(data.emg_hold_avg_f1(tmp_row));
                chord_emg_mat{i}(j,7) = mean(data.emg_hold_avg_f2(tmp_row));
                chord_emg_mat{i}(j,8) = mean(data.emg_hold_avg_f3(tmp_row));
                chord_emg_mat{i}(j,9) = mean(data.emg_hold_avg_f4(tmp_row));
                chord_emg_mat{i}(j,10) = mean(data.emg_hold_avg_f5(tmp_row));
            end
        end
        
        % container for the scaling factors:
        scales = zeros(size(chord_emg_mat{1},2),2);
        
        % looping through sessions:
        for i = 1:size(scales,2)
            % getting single finger patterns of the session:
            tmp = chord_emg_mat{i}(1:10,:);
            
            % calculating the channels norms:
            scales(:,i) = vecnorm(tmp);
        end
        
        varargout{1} = scales;

    case 'inspect_channels_in_natural'
        subject_name = 'subj01';
        sampling_option = 'whole_sampled';
        vararginoptions(varargin,{'subject_name','sampling_option'});

        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_name);

        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        % defining sessions:
        sess = {'sess01','sess02'};
        sess_color = ['k','r'];
        
        % normalizing the natural EMGs:
        for i = 1:length(sess)
            emg_dist{i} = emg_dist{i} ./ scales(:,i)';
        end

        % looping through sessions:
        figure;
        for i = 1:length(sess)
            mean_channels = mean(emg_dist{i},1);
            sem_channels = std(emg_dist{i})/size(emg_dist{i},1);
            
            scatter(1:length(mean_channels),mean_channels,50,sess_color(i),'filled');
            hold on
            errorbar(1:length(mean_channels),mean_channels,sem_channels,'LineWidth',0.5,'color',sess_color(i))
        end
        title('normalized natural EMGs, avg of channels across samples')
        xlim([0,length(mean_channels)+1])
        ylim([0,1])
        xticks(1:length(mean_channels))
        xticklabels(emg_locs_names)

    case 'chord_magnitude_difficulty_model'
        subject_name = 'subj01';
        plot_option = 1;
        vararginoptions(varargin,{'subject_name','plot_option'});

        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};

        % loading subject data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        data = getrow(data,data.sn==str2double(subject_name(end-1:end)));
        
        % calculating avg chord patterns:
        [chord_emg_mat, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name,'plot_option',0,'normalize_channels',1);
        
        % getting avg mean deviation of chords:
        chords_mean_dev = [];
        for j = 1:length(sess)
            for i = 1:length(chords)
                row = data.BN>=sess_blocks{j}(1) & data.BN<=sess_blocks{j}(end) & data.trialCorr==1 & data.chordID==chords(i);
                tmp_mean_dev(i) = mean(data.mean_dev(row));
            end
            chords_mean_dev(:,j) = tmp_mean_dev;
        end
        
        % removing single finger chords:
        % chords_mean_dev(1:10,:) = [];
        % for i = 1:length(sess)
        %     chord_emg_mat{i}(1:10,:) = [];
        % end
        
        % correlations of mean dev and patterns within each session
        corr_sess = zeros(length(sess),1);
        for i = 1:length(sess)
            corr_sess(i) = corr(chords_mean_dev(:,i), mean(chord_emg_mat{i},2));
        end
        
        if (plot_option)
            figure;
            for i = 1:length(sess)
                subplot(2,1,i)
                hold on
                scatter_corr(chords_mean_dev(1:10,i), mean(chord_emg_mat{i}(1:10,:),2), 'r', 'o')
                hold on
                scatter_corr(chords_mean_dev(11:end,i), mean(chord_emg_mat{i}(11:end,:),2), 'k', 'filled')
                title(sprintf('sess %d',i))
                xlabel('avg mean deviation')
                ylabel('avg emg across ch and trial')
                ylim([0,1])
                xlim([0,4])
            end

            figure;
            for i = 1:length(sess)
                subplot(2,1,i)
                hold on
                scatter_corr(chords_mean_dev(:,i), mean(chord_emg_mat{i}(:,:),2), 'k', 'o')
                title(sprintf('sess %d',i))
                xlabel('avg mean deviation')
                ylabel('avg emg across ch and trial')
                ylim([0,1])
                xlim([0,4])
            end
        end
        
        varargout{1} = chords_mean_dev;
        varargout{2} = [mean(chord_emg_mat{1},2),mean(chord_emg_mat{2},2)];

    case 'nSphere_numSamples_vs_radius'
        subject_name = 'subj01';
        d_type = 'Euclidean';
        radius_lim = [0,2];
        n_radius = 100;
        lambda = [];
        sampling_option = 'whole_sampled';
        plot_option = 1;
        vararginoptions(varargin,{'subject_name','d_type','lambda','radius_lim','n_radius','sampling_option','plot_option'})

        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};
        
        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_name);

        % normalizing the natural EMGs:
        for i = 1:length(sess)
            emg_dist{i} = emg_dist{i} ./ scales(:,i)';
        end

        % loading subject data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        data = getrow(data,data.sn==str2double(subject_name(end-1:end)));
        
        % calculating avg chord patterns:
        [chord_emg_mat, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name,'plot_option',0,'normalize_channels',1);       
        
        % vector of radii:
        radius_vec = linspace(radius_lim(1),radius_lim(2),n_radius);

        % container for n samples:
        n_samples = cell(length(sess),1);
        
        % looping through sessions:
        for i = 1:length(sess)
            n_samples{i} = zeros(length(chords),length(radius_vec));
            for j = 1:length(chords)
                fprintf('calculating for session %d , chord %d/%d \n',i,j,length(chords))
                for r = 1:length(radius_vec)
                    [n,~] = nSphere_count_samples(emg_dist{i},chord_emg_mat{i}(j,:),radius_vec(r),'d_type',d_type,'lambda',lambda);
                    n_samples{i}(j,r) = n;
                end
            end
        end

        % plot:
        if plot_option
            figure;
            for i = 1:length(sess)
                subplot(2,2,i);
                % single finger chords
                plot(radius_vec.^size(chord_emg_mat{i},2), n_samples{i}(1:10,:))
                title(sprintf('sess %d , single finger',i))
                xlabel(sprintf('r^{%d}',size(chord_emg_mat{i},2)))
                ylim([0,size(emg_dist{i},1)])
    
                subplot(2,2,i+2);
                % multi finger chords
                plot(radius_vec.^size(chord_emg_mat{i},2), n_samples{i}(11:end,:))
                title(sprintf('sess %d , multi finger',i))
                xlabel(sprintf('r^{%d}',size(chord_emg_mat{i},2)))
                ylim([0,size(emg_dist{i},1)])
            end

            figure;
            for i = 1:length(sess)
                subplot(1,2,i);
                imagesc(n_samples{i})
                title(sprintf('sess %d',i))

                colorbar
                
                % plot settings:
                ax = gca;
                
                set(ax,'YTick',(1:size(chord_emg_mat{i},1)))
                set(ax,'YTickLabel',chords)
                
                set(ax,'XTick', (1:10:length(radius_vec)))
                set(ax,'XTickLabel',round(radius_vec(1:10:end),2))
                
                set(gca,'YDir','reverse')

                xlabel('radius')
            end
        end
        
        varargout{1} = n_samples;
        varargout{2} = radius_vec;

    case 'nSphere_model'
        d_type = 'Euclidean';
        lambda = [];
        n_thresh = 10;
        sampling_option = 'whole_sampled';
        plot_option = 1;
        vararginoptions(varargin,{'d_type','lambda','sampling_option','n_thresh','plot_option'})
        
        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};
        
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        
        % subject numbers:
        subjects = unique(data.sn);
        
        % tranforming subject numbers to subject names:
        subject_names = strcat('subj',num2str(subjects,'%02.f'));

        % container for the dataframe:
        C = [];
        for sn = 1:length(subjects)
            % set natural EMG file name:
            file_name = fullfile(project_path, 'analysis', ['natChord_' subject_names(sn,:) '_emg_natural_' sampling_option '.mat']);
            
            % loading natural EMG dists:
            emg_dist = load(file_name);
            emg_dist = emg_dist.emg_natural_dist;
            
            % scaling factors:
            scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_names(sn,:));
    
            % normalizing the natural EMGs:
            for i = 1:length(sess)
                emg_dist{i} = emg_dist{i} ./ scales(:,i)';
            end

            [avg_patterns, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_names(sn,:),'plot_option',0,'normalize_channels',1);

            % looping through sessions:
            for i = 1:length(sess)
                % container for the each session's dataframe:
                tmp = [];

                % getting avg mean deviation of chords:
                chords_mean_dev = zeros(length(chords),1);
                for j = 1:length(chords)
                    row = data.sn==subjects(sn) & data.BN>=sess_blocks{i}(1) & data.BN<=sess_blocks{i}(end) & data.trialCorr==1 & data.chordID==chords(j);
                    chords_mean_dev(j) = mean(data.mean_dev(row));
                end
    
                for j = 1:length(chords)
                    % sorted distances from the natural dist:
                    d = get_d_from_natural(avg_patterns{i}(j,:)',emg_dist{i},'d_type',d_type,'lambda',lambda);
    
                    % storing the information:
                    tmp.sn(j,1) = subjects(sn);
                    tmp.sess(j,1) = i;
                    tmp.chordID(j,1) = chords(j);
                    tmp.MD(j,1) = chords_mean_dev(j);
                    tmp.thresh(j,1) = n_thresh;
                    tmp.d(j,1) = d(n_thresh);
                    tmp.slope(j,1) = n_thresh/d(n_thresh)^10;
                    tmp.log_slope(j,1) = log(n_thresh/d(n_thresh)^10);
                end

                C = addstruct(C,tmp,'row','force');
            end
        end

        if (plot_option)
            for sn = 1:length(subjects)
                figure;
                for i = 1:length(sess)
                    x = C.MD(C.sn==subjects(sn) & C.sess==i);
                    y = C.slope(C.sn==subjects(sn) & C.sess==i);
                    subplot(1,2,i)
                    hold on
                    scatter_corr(y(1:10) ,x(1:10), 'r', 'o')
                    hold on
                    scatter_corr(y(11:end), x(11:end), 'k', 'o')
                    title(sprintf('Slopes , thresh = %d  , sess %d',C.thresh(1),i))
                    xlabel('slopes')
                    ylabel('avg MD')
                    ylim([0,4])
                end
                sgtitle(sprintf(subject_names(sn,:)))
                
                figure;
                for i = 1:length(sess)
                    x = C.MD(C.sn==subjects(sn) & C.sess==i);
                    y = C.d(C.sn==subjects(sn) & C.sess==i);
                    subplot(1,2,i)
                    hold on
                    scatter_corr(y(1:10), x(1:10), 'r', 'filled')
                    hold on
                    scatter_corr(y(11:end), x(11:end), 'k', 'filled')
                    title(sprintf('thresh = %d , sess %d',C.thresh(1),i))
                    xlabel('d at count threshold')
                    ylabel('avg MD')
                    ylim([0,4])
                end
                sgtitle(subject_names(sn,:))
    
                figure;
                for i = 1:length(sess)
                    x = C.MD(C.sn==subjects(sn) & C.sess==i);
                    y = C.log_slope(C.sn==subjects(sn) & C.sess==i);
                    subplot(1,2,i)
                    hold on
                    scatter_corr(y(1:10,:), x(1:10), 'r', 'filled')
                    hold on
                    scatter_corr(y(11:end,:), x(11:end), 'k', 'filled')
                    title(sprintf('log(slope) , thresh = %d , sess %d',C.thresh(1),i))
                    xlabel('log(slope)')
                    ylabel('avg mean deviation')
                    ylim([0,4])
                end
                sgtitle(subject_names(sn,:))
            end
        end

        varargout{1} = C;

    case 'chord_distance_matrix'
        subject_name = 'subj01';
        normalize_channels = 1;
        vararginoptions(varargin,{'subject_name','normalize_channels'})
        
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        
        % calculating avg chord patterns:
        [chord_emg_mat, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name,'plot_option',0,'normalize_channels',normalize_channels);
        
        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};

        % getting avg mean deviation of chords:
        chords_mean_dev = [];
        c = [];
        for j = 1:length(sess)
            for i = 1:length(chords)
                row = data.BN>=sess_blocks{j}(1) & data.BN<=sess_blocks{j}(end) & data.trialCorr==1 & data.chordID==chords(i) & data.sn==str2double(subject_name(end-1:end));
                tmp_mean_dev(i) = mean(data.mean_dev(row));
            end
            chords_mean_dev(:,j) = tmp_mean_dev;
        end
        
        % emg locations:
        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];
        
        % distance between muscle patterns:
        figure; 
        for i = 1:length(sess)
            d = squareform(pdist(chord_emg_mat{i}));

            subplot(1,2,i)
            imagesc(d)
            colorbar
            
            % plot settings:
            ax = gca;
            
            set(ax,'YTick',(1:size(chords,1)))
            set(ax,'YTickLabel',chords)
            
            set(ax,'XTick', (1:size(chords,1)))
            set(ax,'XTickLabel',chords)
            
            set(gca,'YDir','reverse')
            title(sprintf('emg pattern distance , %s',sess{i}))
        end
        
        figure;
        for i = 1:length(sess)

            % distance of multi finger chords muscle patterns:
            d_patterns = pdist(chord_emg_mat{i}(11:end,:))';

            % difference of mean devs to the power of 2:
            d_MD = pdist(chords_mean_dev(11:end,i))';

            subplot(1,2,i)
            hold on
            scatter(d_MD,d_patterns, 50, 'k', 'filled')
            hold on 
            plot([0,max([max(d_MD),max(d_patterns)])+0.5], [0,max([max(d_MD),max(d_patterns)])+0.5], 'r', 'LineWidth',0.5)
            xlim([0,max([max(d_MD),max(d_patterns)])+0.5])
            ylim([0,max([max(d_MD),max(d_patterns)])+0.5])
            xlabel('difference of MD')
            ylabel('difference of muscle patterns')
            
            title(sess{i})
        end

        % varargout{1} = pdist(chord_emg_mat{1}')

    case 'natural_distance_RDM'
        % handling input arguments:
        sampling_option = 'whole_sampled';
        subject_name = 'subj01';
        vararginoptions(varargin,{'subject_name','sampling_option'});

        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_name);

        sess = {'sess01','sess02'};

        % normalizing the natural EMGs:
        for i = 1:length(sess)
            emg_dist{i} = emg_dist{i} ./ scales(:,i)';
        end
        
        % distance of EMG channels:
        d_emg = [];
        for i = 1:length(sess)
            d_emg{i} = squareform(pdist(emg_dist{i}'));
        end
        
        figure;
        imagesc(d_emg{1})
        colormap('hot')
        colorbar
        figure;
        imagesc(d_emg{2})
        colormap('hot')
        colorbar

    otherwise
        error('The analysis you entered does not exist!')
end



