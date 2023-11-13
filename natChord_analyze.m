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
        vararginoptions(varargin,{'subject_name','smoothing_win_length'});
        
        % if a cell containing multiple subjects was given:
        if (iscell(subject_name))
            for i = 1:length(subject_name)
                efc1_subj(subject_name{i},'smoothing_win_length',smoothing_win_length)
            end
        % if a single subject as a char was given:
        else
            efc1_subj(subject_name,'smoothing_win_length',smoothing_win_length);
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

    case 'avg_chord_patterns'
        % handling input arguments:
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

        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        % looping through sessions:
        figure;
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
            end
            
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

        fprintf("Corr of pattern matrices = %.4f\n",corr2(chord_emg_mat{1},chord_emg_mat{2}))
        varargout{1} = chord_emg_mat;
        varargout{2} = chords;
    
    case 'visualize_natural_emg'
        % handling input arguments:
        sampling_option = 'whole';
        subject_name = 'subj01';
        vararginoptions(varargin,{'subject_name','sampling_option'});

        % defining sessions:
        sess = {'sess01','sess02'};
        sess_blocks = {1:5,6:10};

        % set file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % loading subject data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        
        % calculating avg chord patterns:
        [chord_emg_mat, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name);
        
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
                else
                    scatter3(chord_emg_mat{i}(j,dims(1)), chord_emg_mat{i}(j,dims(2)), chord_emg_mat{i}(j,dims(3)), 100, 'filled', 'MarkerFaceColor', c(j,:))
                end
            end
            legend(num2str(chords(11:end)))

        end
        colorbar;

        figure;
        for i = 1:length(sess)
            subplot(1,2,i)
            
            % pPCA on the natural dist, gets the first 3 dims:
            [COEFF,SCORE,LATENT] = pca(emg_dist{i},'NumComponents',3);

            % scatter 3D natural EMG dist:
            scatter3(SCORE(:,1), SCORE(:,2), SCORE(:,3), 10, 'filled', 'MarkerFaceColor', [0.6,0.6,0.6]);
            xlabel(sprintf('dim 1, var = %.2f',LATENT(1)))
            ylabel(sprintf('dim 2, var = %.2f',LATENT(2)))
            zlabel(sprintf('dim 3, var = %.2f',LATENT(3)))
            title(['pPCA ' sess{i}])
            
            hold all;

            % mapping mean devs to colormap:
            c = map2color(chords_mean_dev(:,i), autumn);

            % put the transformed avg chord patterns on the plot
            for j = 1:size(chord_emg_mat{i},1)
                % in case of single finger chords:
                if (j <= 10)
                    tmp = chord_emg_mat{i}(j,:) * COEFF;
                    scatter3(tmp(1), tmp(2), tmp(3), 100, 'k', 'filled')
                else
                    tmp = chord_emg_mat{i}(j,:) * COEFF;
                    scatter3(tmp(1), tmp(2), tmp(3), 100, 'filled', 'MarkerFaceColor', c(j,:))
                end
            end

        end
        
        
    otherwise
        error('The analysis you entered does not exist!')
end



