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
        files = dir(fullfile(usr_path, 'Desktop', 'Projects', 'EFC1', 'analysis', 'efc1_*_raw.tsv'));
        movFiles = dir(fullfile(usr_path, 'Desktop', 'Projects', 'EFC1', 'analysis', 'efc1_*_mov.mat'));
        
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

        dsave(fullfile(usr_path,'Desktop','Projects','EFC1','analysis','efc1_all.tsv'),ANA);

    case 'chord_avg_emg_visualize'
        % handling input arguments:
        subject_name = 'subj01';
        vararginoptions(varargin,{'subject_name'});

        % loading data:
        data = dload(fullfile(project_path, 'analysis', ['natChord_' subject_name '_raw.tsv']));
        
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

    
    case 'visualize_natural_emg'
        
        sampling_option = 'whole';
        vararginoptions(varargin,{'subject_name','sampling_option'});

        % loading natural EMG dists:

    
    otherwise
        error('The analysis you entered does not exist!')
end



