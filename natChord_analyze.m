function varargout=natChord_analyze(what, varargin)

addpath('functions/')

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);
% project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');
baseDir = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

% colors:
colors_red = [[255, 219, 219] ; [255, 146, 146] ; [255, 73, 73] ; [255, 0, 0] ; [182, 0, 0]]/255;
colors_gray = ['#d3d3d3' ; '#b9b9b9' ; '#868686' ; '#6d6d6d' ; '#535353'];
colors_blue = ['#dbecff' ; '#a8d1ff' ; '#429bff' ; '#0f80ff' ; '#0067db'];
colors_cyan = ['#adecee' ; '#83e2e5' ; '#2ecfd4' ; '#23a8ac' ; '#1b7e81'];
colors_random = ['#773344' ; '#E3B5A4' ; '#83A0A0' ; '#0B0014' ; '#D44D5C'];

colors_blue = hex2rgb(colors_blue);
colors_gray = hex2rgb(colors_gray);
colors_random = hex2rgb(colors_random);
colors_cyan = hex2rgb(colors_cyan);

colors_bold = [colors_red(5,:) ; colors_blue(5,:) ; colors_cyan(5,:) ; colors_gray(5,:) ; colors_random(5,:)];
colors_fade = [colors_red(3,:) ; colors_blue(3,:) ; colors_cyan(3,:) ; colors_gray(3,:) ; colors_random(1,:)];

% figure properties:
my_font.xlabel = 10;
my_font.ylabel = 10;
my_font.title = 11;
my_font.tick_label = 8;
my_font.legend = 8;
my_font.conf_tick_label = 28;
my_font.conf_label = 30;
my_font.conf_legend = 32;
my_font.conf_title = 36;

% conference fig settings:
conf.err_width = 3;
conf.line_width = 8;
conf.marker_size = 350;
conf.horz_line_width = 6;
conf.axis_width = 3;
conf.bar_line_width = 2.5;
conf.bar_width = 2;

switch (what)
    case 'subject_routine'
        % handling input arguments:
        subject_name = 'subj01';
        smoothing_win_length = 30;          % smoothing for force signals
        Fstop_lpf = 40;
        sampling_option = 'whole_sampled';  % sampling option for the natural EMG data
        natural_window_size = 20;          % wn size for natural EMGs sampling in ms.
        wn_spacing = 10;                     % the spacing between windows for the whole_samlped sampling option (i.e. how many windows to skip)
        natural_window_type = 'Rect';       % window shape for the natural EMG sampling. 'Rect' or 'Gaussian'
        vararginoptions(varargin,{'subject_name','smoothing_win_length','Fstop_lpf', ...
                          'sampling_option','natural_window_size','natural_window_type','wn_spacing'});
        
        % if a cell containing multiple subjects was given:
        if (iscell(subject_name))
            for i = 1:length(subject_name)
                natChord_subj(subject_name{i},'smoothing_win_length',smoothing_win_length, 'Fstop_lpf', Fstop_lpf, ...
                             'sampling_option',sampling_option,'natural_window_size',natural_window_size, ...
                             'natural_window_type',natural_window_type,'wn_spacing',wn_spacing);
            end
        % if a single subject as a char was given:
        else
            natChord_subj(subject_name,'smoothing_win_length',smoothing_win_length, 'Fstop_lpf', Fstop_lpf, ...
                         'sampling_option',sampling_option,'natural_window_size',natural_window_size, ...
                         'natural_window_type',natural_window_type,'wn_spacing',wn_spacing);
        end
    
    case 'make_all_dataframe'
        % Calculate RT, MT, MD for each trial of each subejct
        % and create a struct without the mov signals and save it as a
        % single struct called natChord_all.mat
        
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

            diff_force_f1 = zeros(length(tmp_data.BN),1);
            diff_force_f2 = zeros(length(tmp_data.BN),1);
            diff_force_f3 = zeros(length(tmp_data.BN),1);
            diff_force_f4 = zeros(length(tmp_data.BN),1);
            diff_force_f5 = zeros(length(tmp_data.BN),1);

            force_f1 = zeros(length(tmp_data.BN),1);
            force_f2 = zeros(length(tmp_data.BN),1);
            force_f3 = zeros(length(tmp_data.BN),1);
            force_f4 = zeros(length(tmp_data.BN),1);
            force_f5 = zeros(length(tmp_data.BN),1);
            force_e1 = zeros(length(tmp_data.BN),1);
            force_e2 = zeros(length(tmp_data.BN),1);
            force_e3 = zeros(length(tmp_data.BN),1);
            force_e4 = zeros(length(tmp_data.BN),1);
            force_e5 = zeros(length(tmp_data.BN),1);

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
                    % calculate and store rt and mt:
                    [rt_tmp(j),mt_tmp(j),first_finger_tmp(j)] = calculate_rt_mt(tmp_mov{j}, tmp_data.chordID(j), ...
                                                                tmp_data.baselineTopThresh(j), tmp_data.RT(j), ...
                                                                tmp_data.fGain1(j), tmp_data.fGain2(j), tmp_data.fGain3(j), ...
                                                                tmp_data.fGain4(j), tmp_data.fGain5(j));

                    % average force:
                    idx_completion = find(tmp_mov{j}(:,1)==3);
                    idx_completion = idx_completion(end);
                    diff_force_f1(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,14));
                    diff_force_f2(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,15));
                    diff_force_f3(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,16));
                    diff_force_f4(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,17));
                    diff_force_f5(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,18));

                    force_f1(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,9));
                    force_f2(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,10));
                    force_f3(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,11));
                    force_f4(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,12));
                    force_f5(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,13));
                    force_e1(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,4));
                    force_e2(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,5));
                    force_e3(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,6));
                    force_e4(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,7));
                    force_e5(j) = mean(tmp_mov{j}(idx_completion-299:idx_completion,8));
                
                % if trial was incorrect:
                else
                    diff_force_f1(j) = -1;
                    diff_force_f2(j) = -1;
                    diff_force_f3(j) = -1;
                    diff_force_f4(j) = -1;
                    diff_force_f5(j) = -1;

                    force_f1(j) = -1;
                    force_f2(j) = -1;
                    force_f3(j) = -1;
                    force_f4(j) = -1;
                    force_f5(j) = -1;
                    force_e1(j) = -1;
                    force_e2(j) = -1;
                    force_e3(j) = -1;
                    force_e4(j) = -1;
                    force_e5(j) = -1;

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
            tmp_data.MD = mean_dev_tmp;

            sess = (tmp_data.BN<=10) + 2*(tmp_data.BN>=11);
            tmp_data.sess = sess;

            tmp_data.diff_force_f1 = diff_force_f1;
            tmp_data.diff_force_f2 = diff_force_f2;
            tmp_data.diff_force_f3 = diff_force_f3;
            tmp_data.diff_force_f4 = diff_force_f4;
            tmp_data.diff_force_f5 = diff_force_f5;

            tmp_data.force_f1 = force_f1;
            tmp_data.force_f2 = force_f2;
            tmp_data.force_f3 = force_f3;
            tmp_data.force_f4 = force_f4;
            tmp_data.force_f5 = force_f5;
            tmp_data.force_e1 = force_e1;
            tmp_data.force_e2 = force_e2;
            tmp_data.force_e3 = force_e3;
            tmp_data.force_e4 = force_e4;
            tmp_data.force_e5 = force_e5;
            
            % adding subject data to ANA:
            ANA=addstruct(ANA,tmp_data,'row','force');
        end
        % adding number of active fingers:
        ANA.num_fingers = get_num_active_fingers(ANA.chordID);

        dsave(fullfile(project_path,'analysis','natChord_all.tsv'),ANA);

    case 'make_chord_dataframe'
        % fields:
        % sn, sess, chordID, num_trials, num_fingers, MD, MT, RT, MD_std, MT_std, RT_std  
        
        % load trial dataframe:
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        subjects = unique(data.sn);
        sess = data.sess;
        chords = unique(data.chordID);
        % sorting chords in an arbitrary way:
        chords_sorted = [19999, 91999, 99199, 99919, 99991, 29999, 92999, 99299, 99929, 99992];
        [ind,~] = find(chords == chords_sorted);
        chords = [chords_sorted'; chords(setdiff(1:length(chords),ind))];
        
        n = get_num_active_fingers(chords);

        % container to hold the dataframe:
        ANA = [];
        % loop on subjects:
        for i = 1:length(subjects)
            subj_sess = unique(data.sess(data.sn==subjects(i)));
            subj_sess = [1]; % Only analyze the data for session 1
            % loop on sess:
            for j = 1:length(subj_sess)
                % loop on chords:
                for k = 1:length(chords)
                    tmp = [];
                    
                    tmp.sn = subjects(i);
                    tmp.sess = j;
                    tmp.chordID = chords(k);
                    
                    row = data.sn==subjects(i) & sess==j & data.chordID==chords(k) & data.trialCorr==1;
                    tmp.num_trials = sum(row);
                    tmp.num_fingers = n(k);
                    tmp.MD = mean(data.MD(row));
                    tmp.MT = mean(data.MT(row));
                    tmp.RT = mean(data.RT(row));
                    tmp.MD_std = std(data.MD(row));
                    tmp.MT_std = std(data.MT(row));
                    tmp.RT_std = std(data.RT(row));

                    tmp.diff_force_f1 = mean(data.diff_force_f1(row));
                    tmp.diff_force_f2 = mean(data.diff_force_f2(row));
                    tmp.diff_force_f3 = mean(data.diff_force_f3(row));
                    tmp.diff_force_f4 = mean(data.diff_force_f4(row));
                    tmp.diff_force_f5 = mean(data.diff_force_f5(row));

                    tmp.force_f1 = mean(data.force_f1(row));
                    tmp.force_f2 = mean(data.force_f2(row));
                    tmp.force_f3 = mean(data.force_f3(row));
                    tmp.force_f4 = mean(data.force_f4(row));
                    tmp.force_f5 = mean(data.force_f5(row));
                    tmp.force_e1 = mean(data.force_e1(row));
                    tmp.force_e2 = mean(data.force_e2(row));
                    tmp.force_e3 = mean(data.force_e3(row));
                    tmp.force_e4 = mean(data.force_e4(row));
                    tmp.force_e5 = mean(data.force_e5(row));
                    
                    tmp.emg_hold_avg_e1 = mean(data.emg_hold_avg_e1(row));
                    tmp.emg_hold_avg_e2 = mean(data.emg_hold_avg_e2(row));
                    tmp.emg_hold_avg_e3 = mean(data.emg_hold_avg_e3(row));
                    tmp.emg_hold_avg_e4 = mean(data.emg_hold_avg_e4(row));
                    tmp.emg_hold_avg_e5 = mean(data.emg_hold_avg_e5(row));
                    tmp.emg_hold_avg_f1 = mean(data.emg_hold_avg_f1(row));
                    tmp.emg_hold_avg_f2 = mean(data.emg_hold_avg_f2(row));
                    tmp.emg_hold_avg_f3 = mean(data.emg_hold_avg_f3(row));
                    tmp.emg_hold_avg_f4 = mean(data.emg_hold_avg_f4(row));
                    tmp.emg_hold_avg_f5 = mean(data.emg_hold_avg_f5(row));

                    ANA = addstruct(ANA,tmp,'row','force');
                end
            end
        end
        ANA.scale_e1 = zeros(length(ANA.chordID),1);
        ANA.scale_e2 = zeros(length(ANA.chordID),1);
        ANA.scale_e3 = zeros(length(ANA.chordID),1);
        ANA.scale_e4 = zeros(length(ANA.chordID),1);
        ANA.scale_e5 = zeros(length(ANA.chordID),1);
        ANA.scale_f1 = zeros(length(ANA.chordID),1);
        ANA.scale_f2 = zeros(length(ANA.chordID),1);
        ANA.scale_f3 = zeros(length(ANA.chordID),1);
        ANA.scale_f4 = zeros(length(ANA.chordID),1);
        ANA.scale_f5 = zeros(length(ANA.chordID),1);
        
        % chords = [19999,91999,99199,99919,99991,29999,92999,99299,99929,99992];
        % adding EMG scale factors to dataframe:
        % loop on subjects:
        for i = 1:length(subjects)
            % loop on sess:
            for j = 1:length(unique(sess))
                row = ANA.sn == subjects(i) & ANA.sess == j;
                e1_profile = [];
                e2_profile = [];
                e3_profile = [];
                e4_profile = [];
                e5_profile = [];
                f1_profile = [];
                f2_profile = [];
                f3_profile = [];
                f4_profile = [];
                f5_profile = [];
                % activity profile of each channel across chords:
                for k = 1:length(chords)
                    e1_profile = [e1_profile ; ANA.emg_hold_avg_e1(row & ANA.chordID==chords(k))];
                    e2_profile = [e2_profile ; ANA.emg_hold_avg_e2(row & ANA.chordID==chords(k))];
                    e3_profile = [e3_profile ; ANA.emg_hold_avg_e3(row & ANA.chordID==chords(k))];
                    e4_profile = [e4_profile ; ANA.emg_hold_avg_e4(row & ANA.chordID==chords(k))];
                    e5_profile = [e5_profile ; ANA.emg_hold_avg_e5(row & ANA.chordID==chords(k))];
                    f1_profile = [f1_profile ; ANA.emg_hold_avg_f1(row & ANA.chordID==chords(k))];
                    f2_profile = [f2_profile ; ANA.emg_hold_avg_f2(row & ANA.chordID==chords(k))];
                    f3_profile = [f3_profile ; ANA.emg_hold_avg_f3(row & ANA.chordID==chords(k))];
                    f4_profile = [f4_profile ; ANA.emg_hold_avg_f4(row & ANA.chordID==chords(k))];
                    f5_profile = [f5_profile ; ANA.emg_hold_avg_f5(row & ANA.chordID==chords(k))];
                end
                % norm of each channel is the scale:
                ANA.scale_e1(row) = norm(e1_profile);
                ANA.scale_e2(row) = norm(e2_profile);
                ANA.scale_e3(row) = norm(e3_profile);
                ANA.scale_e4(row) = norm(e4_profile);
                ANA.scale_e5(row) = norm(e5_profile);
                ANA.scale_f1(row) = norm(f1_profile);
                ANA.scale_f2(row) = norm(f2_profile);
                ANA.scale_f3(row) = norm(f3_profile);
                ANA.scale_f4(row) = norm(f4_profile);
                ANA.scale_f5(row) = norm(f5_profile);
            end
        end

        % Adding measures from EFC1 dataset: 
        ANA.MD_efc = zeros(size(ANA.MD));
        ANA.MT_efc = zeros(size(ANA.MD));
        ANA.RT_efc = zeros(size(ANA.MD));
        subject_mapping = [subjects , [1,5,2,4,7,15,17,11,8,16]'];
        data_efc1 = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        chords = unique(ANA.chordID);
        for i = 1:length(subjects)
            for j = 1:length(chords)
                % select part of the data from efc1:
                subject_efc = subject_mapping(i,2);
                tmp = getrow(data_efc1,data_efc1.sn==subject_efc & data_efc1.sess>=3 & data_efc1.chordID==chords(j));
                
                ANA.MD_efc(ANA.sn==subjects(i) & ANA.chordID==chords(j)) = mean(tmp.MD,'omitmissing');
                ANA.MT_efc(ANA.sn==subjects(i) & ANA.chordID==chords(j)) = mean(tmp.MT,'omitmissing');
                ANA.RT_efc(ANA.sn==subjects(i) & ANA.chordID==chords(j)) = mean(tmp.RT,'omitmissing');
            end
        end
        
        dsave(fullfile(project_path,'analysis','natChord_chord.tsv'),ANA);

    case 'make_analysis_dataframe'
        % load chord dataframe:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        
        % nSphere model:
        [C,~] = natChord_analyze('nSphere_model','d_type','project_to_nSphere','sampling_option','whole_thresholded','plot_option',0);

        data.log_slope = C.log_slope;
        data.log_slope_n = C.log_slope_n;
        data.d = C.d;
        
        % magnitude model:
        [~,out] = natChord_analyze('chord_magnitude','plot_option',0);
        data.magnitude = out.mag;
        data.magnitude_n = out.mag_n;
            
        % coactivation model:
        [~,out] = natChord_analyze('chord_coact','plot_option',0);
        data.coact = out.coact;

        % emg channel values:
        C = natChord_analyze('EMG_prewhitening_matrix','plot_option',0);
        subj = unique(C.sn);
        emg_f1 = [];
        emg_f2 = [];
        emg_f3 = [];
        emg_f4 = [];
        emg_f5 = [];
        emg_e1 = [];
        emg_e2 = [];
        emg_e3 = [];
        emg_e4 = [];
        emg_e5 = [];
        emg_f1_pw = [];
        emg_f2_pw = [];
        emg_f3_pw = [];
        emg_f4_pw = [];
        emg_f5_pw = [];
        emg_e1_pw = [];
        emg_e2_pw = [];
        emg_e3_pw = [];
        emg_e4_pw = [];
        emg_e5_pw = [];
        for i = 1:length(subj)
            sess = unique(C.sess(C.sn==subj(i)));
            for j = 1:length(sess)
                emg_patterns = C.pattern{C.sn==subj(i) & C.sess==sess(j)};
                emg_patterns_pw = C.pattern_prewhitened{C.sn==subj(i) & C.sess==sess(j)};
                
                emg_f1 = [emg_f1 ; emg_patterns(:,6)];
                emg_f2 = [emg_f2 ; emg_patterns(:,7)];
                emg_f3 = [emg_f3 ; emg_patterns(:,8)];
                emg_f4 = [emg_f4 ; emg_patterns(:,9)];
                emg_f5 = [emg_f5 ; emg_patterns(:,10)];
                emg_e1 = [emg_e1 ; emg_patterns(:,1)];
                emg_e2 = [emg_e2 ; emg_patterns(:,2)];
                emg_e3 = [emg_e3 ; emg_patterns(:,3)];
                emg_e4 = [emg_e4 ; emg_patterns(:,4)];
                emg_e5 = [emg_e5 ; emg_patterns(:,5)];
                emg_f1_pw = [emg_f1_pw ; emg_patterns_pw(:,6)];
                emg_f2_pw = [emg_f2_pw ; emg_patterns_pw(:,7)];
                emg_f3_pw = [emg_f3_pw ; emg_patterns_pw(:,8)];
                emg_f4_pw = [emg_f4_pw ; emg_patterns_pw(:,9)];
                emg_f5_pw = [emg_f5_pw ; emg_patterns_pw(:,10)];
                emg_e1_pw = [emg_e1_pw ; emg_patterns_pw(:,1)];
                emg_e2_pw = [emg_e2_pw ; emg_patterns_pw(:,2)];
                emg_e3_pw = [emg_e3_pw ; emg_patterns_pw(:,3)];
                emg_e4_pw = [emg_e4_pw ; emg_patterns_pw(:,4)];
                emg_e5_pw = [emg_e5_pw ; emg_patterns_pw(:,5)];
            end
        end
        data.emg_f1 = emg_f1;
        data.emg_f2 = emg_f2;
        data.emg_f3 = emg_f3;
        data.emg_f4 = emg_f4;
        data.emg_f5 = emg_f5;
        data.emg_e1 = emg_e1;
        data.emg_e2 = emg_e2;
        data.emg_e3 = emg_e3;
        data.emg_e4 = emg_e4;
        data.emg_e5 = emg_e5;

        data.emg_f1_pw = emg_f1_pw;
        data.emg_f2_pw = emg_f2_pw;
        data.emg_f3_pw = emg_f3_pw;
        data.emg_f4_pw = emg_f4_pw;
        data.emg_f5_pw = emg_f5_pw;
        data.emg_e1_pw = emg_e1_pw;
        data.emg_e2_pw = emg_e2_pw;
        data.emg_e3_pw = emg_e3_pw;
        data.emg_e4_pw = emg_e4_pw;
        data.emg_e5_pw = emg_e5_pw;
        
        dsave(fullfile(project_path,'analysis','natChord_analysis.tsv'),data);

    case 'make_natural_dist'
        subject_name = 'subj01';
        fs_emg = 2148.1481;                     % EMG sampling rate in Hz
        Fstop_lpf = 40;
        natural_window_size = 20;       % window size to sample natural EMG
        sampling_option = 'whole_thresholded';      % sampling option to select windows from natural EMGs.
        natural_window_type = 'Rect';   % sampling window type for natural EMGs.
        wn_spacing = 10;                % sampling spacing for the 'whole_sampled' option.
        vararginoptions(varargin,{'subject_name','Fstop_lpf', ...
                                  'sampling_option','natural_window_size','natural_window_type','wn_spacing'});

        data = dload(fullfile(project_path,'analysis','natChord_all.tsv'));
        
        % EMG filters:
        % Define filter parameters
        % Design the Butterworth filter
        [b, a] = butter(6, Fstop_lpf/(fs_emg/2), 'low'); % 6th order Butterworth filter
        % Convert to second-order sections (SOS) format
        [sos, g] = tf2sos(b, a);
                
        % if a cell containing multiple subjects was given:
        if (iscell(subject_name))
            for i = 1:length(subject_name)
                sess = unique(data.sess(data.sn==str2double(subject_name{i}(end-1:end))));
                sess = [1];
                sess_cell = cellfun(@(x) ['sess', sprintf('%02d', x)], num2cell(sess), 'UniformOutput', false);

                % Preprocessing and dealing with the natural EMGs:
                fprintf("Processing natural EMG data...\n\n")
                make_natural_emg(subject_name{i}, sess_cell, fs_emg, sos, g,natural_window_type,natural_window_size,sampling_option,wn_spacing);
            end
        % if a single subject as a char was given:
        else
            sess = unique(data.sess(data.sn==str2double(subject_name(end-1:end))));
            sess = [1];
            sess_cell = cellfun(@(x) ['sess', sprintf('%02d', x)], num2cell(sess), 'UniformOutput', false);

            % Preprocessing and dealing with the natural EMGs:
            fprintf("Processing natural EMG data...\n\n")
            make_natural_emg(subject_name, sess_cell, fs_emg, sos, g,natural_window_type,natural_window_size,sampling_option,wn_spacing);
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
        plot_option = 1;            % flag to plot the avg chord patterns.
        normalize_channels = 1;     % flag to whether normalize the channels by their norms or not.
        sess = 1;
        vararginoptions(varargin,{'subject_name','plot_option','normalize_channels','sess'});
        

        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        data = getrow(data, data.sess == sess & data.sn == str2double(subject_name(end-1:end)));
        chords = unique(data.chordID);
        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];
        
        % emg scaling factors:
        scales = [data.scale_e1(1),data.scale_e2(1),data.scale_e3(1),data.scale_e4(1),data.scale_e5(1), ...
                  data.scale_f1(1),data.scale_f2(1),data.scale_f3(1),data.scale_f4(1),data.scale_f5(1)];

        % sorting chords in an arbitrary way:
        chords_sorted = [19999, 91999, 99199, 99919, 99991, 29999, 92999, 99299, 99929, 99992];
        [ind,~] = find(chords == chords_sorted);
        chords = [chords_sorted'; chords(setdiff(1:length(chords),ind))];

        pattern = zeros(length(chords),10);
        % looping through chords of the sessions:
        for i = 1:length(chords)
            tmp(1,1) = [data.emg_hold_avg_e1(data.chordID==chords(i))];
            tmp(1,2) = [data.emg_hold_avg_e2(data.chordID==chords(i))];
            tmp(1,3) = [data.emg_hold_avg_e3(data.chordID==chords(i))];
            tmp(1,4) = [data.emg_hold_avg_e4(data.chordID==chords(i))];
            tmp(1,5) = [data.emg_hold_avg_e5(data.chordID==chords(i))];
            tmp(1,6) = [data.emg_hold_avg_f1(data.chordID==chords(i))];
            tmp(1,7) = [data.emg_hold_avg_f2(data.chordID==chords(i))];
            tmp(1,8) = [data.emg_hold_avg_f3(data.chordID==chords(i))];
            tmp(1,9) = [data.emg_hold_avg_f4(data.chordID==chords(i))];
            tmp(1,10) = [data.emg_hold_avg_f5(data.chordID==chords(i))];
            
            pattern(i,:) = tmp;

            if (normalize_channels)
                pattern(i,:) = pattern(i,:) ./ scales;
            end
        end
        

        if (plot_option)
            figure('Position',[10 10 450 1500]);
            % plotting the chord EMG patterns:
            pcolor([[pattern, zeros(size(pattern,1),1)] ; zeros(1,size(pattern,2)+1)])
            colorbar
            % clim([0 1.5])
            
            % plot settings:
            ax = gca;
            
            set(ax,'YTick',(1:size(pattern,1))+0.5)
            set(ax,'YTickLabel',chords)
            
            set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
            set(ax,'XTickLabel',emg_locs_names)
            
            set(gca,'YDir','reverse')
            title(['sess ' sess])
        end
        
        varargout{1} = pattern;
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

    case 'behavior_trends'
        measure = 'MD';
        conference_fig = 0;
        vararginoptions(varargin,{'measure','conference_fig'})

        % loading data:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));

        % getting the values of measure:
        values = eval(['data.' measure]);
        
        cond_vec = data.num_fingers;
        [sem_subj, ~, ~] = get_sem(values, data.sn, data.sess, cond_vec);

        if ~conference_fig
            % avg trend acorss sessions:
            fig = figure('Position', [500 500 100 200]);
            fontsize(fig, my_font.tick_label, 'points')
    
            errorbar(sem_subj.partitions(sem_subj.cond==1),sem_subj.y(sem_subj.cond==1),sem_subj.sem(sem_subj.cond==1),'LineStyle','none','Color',colors_blue(2,:)); hold on;
            lineplot(data.sess(data.num_fingers==1),values(data.num_fingers==1),'markertype','o','markersize',5,'markerfill',colors_blue(2,:),'markercolor',colors_blue(2,:),'linecolor',colors_blue(2,:),'linewidth',2,'errorbars','');hold on;
            
            errorbar(sem_subj.partitions(sem_subj.cond==3),sem_subj.y(sem_subj.cond==3),sem_subj.sem(sem_subj.cond==3),'LineStyle','none','Color',colors_blue(3,:))
            lineplot(data.sess(data.num_fingers==3),values(data.num_fingers==3),'markertype','o','markersize',5,'markerfill',colors_blue(3,:),'markercolor',colors_blue(3,:),'linecolor',colors_blue(3,:),'linewidth',2,'errorbars','');
           
            errorbar(sem_subj.partitions(sem_subj.cond==5),sem_subj.y(sem_subj.cond==5),sem_subj.sem(sem_subj.cond==5),'LineStyle','none','Color',colors_blue(5,:))
            lineplot(data.sess(data.num_fingers==5),values(data.num_fingers==5),'markertype','o','markersize',5,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',2,'errorbars','');
            
            % legend('single finger','chord','');
            % legend boxoff
            xlabel('sess','FontSize',my_font.xlabel)
            ylabel('','FontSize',my_font.title)
            % ylim([0.2 2.7])
            % ylim([0 2500])
            ylim([0 500])
            xlim([0.8 2.2])
            % h = gca;
            % h.YTick = linspace(h.YTick(1),h.YTick(end),5);
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
        else
            fig = figure('Units','centimeters', 'Position',[15 15 10 20]);
            fontsize(fig, my_font.conf_tick_label, 'points')
            
            errorbar(sem_subj.partitions(sem_subj.cond==1),sem_subj.y(sem_subj.cond==1),sem_subj.sem(sem_subj.cond==1),'LineStyle','none','Color',colors_blue(2,:),'LineWidth',conf.err_width); hold on;
            lineplot(data.sess(data.num_fingers==1),values(data.num_fingers==1),'markertype','o','markersize',12,'markerfill',colors_blue(2,:),'markercolor',colors_blue(2,:),'linecolor',colors_blue(2,:),'linewidth',6,'errorbars','');hold on;
            
            errorbar(sem_subj.partitions(sem_subj.cond==3),sem_subj.y(sem_subj.cond==3),sem_subj.sem(sem_subj.cond==3),'LineStyle','none','Color',colors_blue(3,:),'LineWidth',conf.err_width)
            lineplot(data.sess(data.num_fingers==3),values(data.num_fingers==3),'markertype','o','markersize',12,'markerfill',colors_blue(3,:),'markercolor',colors_blue(3,:),'linecolor',colors_blue(3,:),'linewidth',6,'errorbars','');
           
            errorbar(sem_subj.partitions(sem_subj.cond==5),sem_subj.y(sem_subj.cond==5),sem_subj.sem(sem_subj.cond==5),'LineStyle','none','Color',colors_blue(5,:),'LineWidth',conf.err_width)
            lineplot(data.sess(data.num_fingers==5),values(data.num_fingers==5),'markertype','o','markersize',12,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',6,'errorbars','');

            xlim([0.8 2.2])

            xlabel('days','FontSize',my_font.conf_label)
            ylabel('','FontSize',my_font.conf_label)

            if measure=='MD'
                ylim([0.2 2.7])
            elseif measure=='RT'
                 ylim([150 450])
            elseif measure=='MT'
                ylim([0 3500])
            end
            h = gca;
            h.YTick = linspace(h.YTick(1),h.YTick(end),3);
            h.XAxis.FontSize = my_font.conf_tick_label;
            h.YAxis.FontSize = my_font.conf_tick_label;
            h.LineWidth = conf.axis_width;
            fontname("Arial")
        end

    
    case 'behavior_reliability'
        % getting subject numbers:
        data_nat = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        data_behav = dload(fullfile(project_path, 'analysis', 'efc1_chord.tsv'));

        % selecting the corresponding portion of the data from behavior
        % experiment:
        chords = unique(data_nat.chordID);
        row01 = arrayfun(@(x) ~isempty(intersect(x,chords)), data_behav.chordID);
        row02 = arrayfun(@(x) ~isempty(intersect(x,[1;5])), data_behav.sn);
        row03 = data_behav.sess >= 3;
        data_behav = getrow(data_behav,row01 & row02 & row03);

        % merging data:
        data_behav.sn(data_behav.sn==5) = 2;
        data_behav.sess = data_behav.sess-2;
        data = addstruct(data_nat,data_behav,'row','force');
        data.set = kron([1;2],ones(length(data_nat.sn),1));

        % struct to hold the results:
        C_w = [];

        % within subject:
        subjects = unique(data.sn);
        cond_vec =  double(data.num_fingers>1) + 1;
        conditions = unique(cond_vec);
        cnt = 1;
        for j = 2:length(conditions)
            for i = 1:length(subjects)
                % natural dataset within subject alpha:
                row = data.set==1 & data.sn==subjects(i) & cond_vec==conditions(j);
                X1 = reshape(data.MD(row), [length(data.sess(row))/length(unique(data.sess(row))), length(unique(data.sess(row)))]);
                alpha_nat = cronbach(X1);

                % natural dataset within subject alpha:
                row = data.set==2 & data.sn==subjects(i) & cond_vec==conditions(j);
                X2 = reshape(data.MD(row), [length(data.sess(row))/length(unique(data.sess(row))), length(unique(data.sess(row)))]);
                alpha_behav = cronbach(X2);
                
                % soring chords in datasets to make the two sets
                % correspondant:
                if j>1
                    [~,i_nat] = sort(data.chordID(data.set==2 & data.sess==1 & data.sn==subjects(i) & cond_vec==conditions(j)));
                    X2 = X2(i_nat,:);
                end

                % between dataset within subject:
                alpha = 0;
                for k1 = 1:size(X1,2)
                    for k2 = 1:size(X2,2)
                        alpha = alpha + cronbach([X1(:,k1) , X2(:,k2)])/size(X1,2)/size(X2,2);
                    end
                end
                
                C_w.cond(cnt,1) = conditions(j);
                C_w.sn(cnt,1) = subjects(i);
                C_w.r_nat(cnt,1) = alpha_nat;
                C_w.r_beh(cnt,1) = alpha_behav;
                C_w.r(cnt,1) = alpha;
                cnt = cnt+1;
            end
        end

        % plot:
        fig = figure('Position',[500 500 250 250]);
        fontsize(fig,my_font.tick_label,'points');
        scatter(C_w.sn, C_w.r, 50, 'MarkerEdgeColor', colors_blue(5,:), 'MarkerFaceColor', colors_blue(5,:)); hold on;
        scatter(C_w.sn, sqrt(C_w.r_beh .* C_w.r_nat), 50, 'MarkerEdgeColor', colors_red(5,:), 'MarkerFaceColor', colors_red(5,:))
        xlabel('subject','FontSize',my_font.xlabel)
        ylabel('Cronbach Alpha','FontSize',my_font.ylabel)
        title('Within Subject Reliability')
        xlim([0.7 2.3])
        ylim([0 1.3])
        drawline(1,'dir','horz','color',[0.7 0.7 0.7]); 
        legend('across dataset','within dataset','')
        legend boxoff
        h = gca;
        h.YTick = 0:0.2:1;
        h.XTick = 1:length(subjects);

        varargout{1} = C_w;
        

    case 'behavior_var_decomp'
        measure = 'MD';
        centered = 1;
        cond = 1;
        vararginoptions(varargin,{'measure','centered','cond'})

        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));

        % getting the values of measure:
        values = eval(['data.' measure]);

        if cond == 1
            cond_vec = double(data.num_fingers>1) + 1;
            [v_g, v_gs, v_gse] = reliability_var(values, data.sn, data.sess, 'cond_vec', cond_vec, 'centered', centered);
        else
            [v_g, v_gs, v_gse] = reliability_var(values, data.sn, data.sess, 'centered', centered);
        end
        
        % plot:
        if cond == 1
            y = [];
            fig = figure();
            fontsize(fig,my_font.tick_label,"points")
            for i = 1:length(unique(cond_vec))
                y(i,:) = [v_g{i}/v_gse{i} (v_gs{i}-v_g{i})/v_gse{i} (v_gse{i}-v_gs{i})/v_gse{i}];
                b = bar(i,y(i,:),'stacked','FaceColor','flat');
                b(1).CData = [238, 146, 106]/255;   % global var
                b(2).CData = [36, 168, 255]/255;  % subj var
                b(3).CData = [0.8 0.8 0.8];  % noise var
                b(1).EdgeColor = [1 1 1];
                b(2).EdgeColor = [1 1 1];
                b(3).EdgeColor = [1 1 1];
                b(1).LineWidth = 3;
                b(2).LineWidth = 3;
                b(3).LineWidth = 3;
                hold on
            end
            drawline(1,'dir','horz','color',[0.7 0.7 0.7])
            drawline(0,'dir','horz','color',[0.7 0.7 0.7])
            box off;
            title(['var decomp ' replace(measure,'_',' ')],'FontSize',my_font.title)
            xlabel('num fingers','FontSize',my_font.xlabel)
            ylabel('percent variance','FontSize',my_font.ylabel)
            title(['var decomp ' replace(measure,'_',' ')],'FontSize',my_font.title);
            legend('g','s','e')
            legend boxoff
            ylim([-0.3,1.2])
            h = gca;
            h.XTick = 1:length(unique(cond_vec));
            h.XTickLabel = {'single finger', 'chord'};
        else
            fig = figure();
            fontsize(fig, my_font.tick_label, "points")
            bar(1,1,'FaceColor','flat','EdgeColor',[1,1,1],'LineWidth',4,'CData',[0.8 0.8 0.8])
            hold on
            bar(1,1-(v_gs-v_g)/v_gse,'FaceColor','flat','EdgeColor',[1,1,1],'LineWidth',4,'CData',[36, 168, 255]/255)
            bar(1,v_g/v_gse,'FaceColor','flat','EdgeColor',[1,1,1],'LineWidth',4,'CData',[238, 146, 106]/255)
            drawline(1,'dir','horz','color',[0.8 0.8 0.8])
            legend('e','s','g')
            legend boxoff
            box off
            ylim([0 1.2])
            xticklabels('all chords')
            ylabel('percent variance','FontSize',my_font.ylabel)
            title(['var decomp ' replace(measure,'_',' ')],'FontSize',my_font.title);
        end

        varargout{1} = [v_g, v_gs, v_gse];



    case 'visualize_natural_emg'
        % handling input arguments:
        sampling_option = 'whole_thresholded';
        subject_name = 'subj01';
        normalize_channels = 1;             % flag to whether normalize the channels by their norms or not.
        dimensions = [];                    % dimensions of the natural data to show. by default random dimensions are selected.
        sess = 1;
        ica = 0;
        n_components = 5;
        vararginoptions(varargin,{'subject_name','sampling_option','normalize_channels','dimensions','sess','ica','n_components'});

        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = get_emg_scales(str2double(subject_name(end-1:end)),sess);
        
        if (normalize_channels)
            % normalizing the natural EMGs:
            for i = 1:length(emg_dist.dist)
                emg_dist.dist{i} = emg_dist.dist{i} ./ scales;
            end
        end

        if ica
            for i = 1:length(emg_dist.dist)
                tmp = transform_dist(emg_dist.dist{i},n_components);
                emg_dist.dist{i} = tmp.transformed_dist{1};
                emg_dist.transform_mat(i,1) = tmp.transform_mat;
                emg_dist.mdl(i,1) = tmp.mdl;
            end
        end
        
        % loading subject data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        
        % calculating avg chord patterns:
        [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_name,'plot_option',0,'normalize_channels',normalize_channels);

        % getting avg mean deviation of chords:
        chords_mean_dev = zeros(length(chords),1);
        for i = 1:length(chords)
            row = data.sess==sess & data.chordID==chords(i) & data.sn==str2double(subject_name(end-1:end));
            chords_mean_dev(i) = data.MD(row);
        end
        
        % emg locations:
        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        if ica
            emg_locs_names = 1:n_components;
        end
        
        % select 3 random dimensions:
        dims = randperm(size(pattern,2));
        dims = dims(1:3);

        if ica
            dims = randperm(n_components);
            dims = dims(1:3);
        end
        
        % if user input dimensions:
        if (~isempty(dimensions))
            if length(dimensions) ~= 3
                warning('Input dimensions length must be 3 -> Changed to random dimensions...')
            else
                dims = dimensions;
            end
        end

        % PLOTS - Natural Stats + chord EMG patterns:
        % loop on partitions:
        for part = 1:length(unique(emg_dist.partition))
            figure;

            tmp_dist = emg_dist.dist(emg_dist.sess==sess);
            tmp_dist = tmp_dist{part};

            % scatter 3D natural EMG dist:
            scatter3(tmp_dist(:,dims(1)), tmp_dist(:,dims(2)), tmp_dist(:,dims(3)), 5, 'filled', 'MarkerFaceColor', [0.6,0.6,0.6], 'HandleVisibility','off');
            xlabel(emg_locs_names(dims(1)),'FontSize',my_font.xlabel)
            ylabel(emg_locs_names(dims(2)),'FontSize',my_font.ylabel)
            zlabel(emg_locs_names(dims(3)),'FontSize',my_font.xlabel)
            title([num2str(sess) ' , partition ' num2str(part)],'FontSize',my_font.title)
            hold on;
            
            % mapping mean devs to colormap:
            c = map2color(chords_mean_dev, autumn);
            
            for j = 1:size(pattern,1)
                scatter3(pattern(j,dims(1)), pattern(j,dims(2)), pattern(j,dims(3)), 50, 'filled', 'MarkerFaceColor', c(j,:))
            end

            % legend(num2str(chords(11:end)))
            colorbar;
        end


    case 'inspect_channels_in_natural'
        subject_name = 'subj01';
        sampling_option = 'whole_sampled';
        sess = 1;
        vararginoptions(varargin,{'subject_name','sampling_option','sess'});

        % set natural EMG file name:
        file_name = fullfile(project_path, 'analysis', ['natChord_' subject_name '_emg_natural_' sampling_option '.mat']);
        
        % loading natural EMG dists:
        emg_dist = load(file_name);
        emg_dist = emg_dist.emg_natural_dist;

        % scaling factors:
        scales = get_emg_scales(str2double(subject_name(end-1:end)),sess);

        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];
        
        % normalizing the natural EMGs:
        for i = 1:length(emg_dist.dist)
            emg_dist.dist{i} = emg_dist.dist{i} ./ scales;
        end

        mean_channels = zeros(length(emg_dist.dist),size(emg_dist.dist{1},2));
        for i = 1:length(emg_dist.dist)
            mean_channels(i,:) = mean(emg_dist.dist{i},1);
        end
        % plot:
        fig = figure('Position',[500 500 300 200]);
        fontsize(fig,my_font.tick_label,'points');
        for i = 1:length(unique(emg_dist.sess))
            scatter(1:size(mean_channels,2), mean(mean_channels(emg_dist.sess==i,:),1), 50, 'MarkerEdgeColor','none','MarkerFaceColor',colors_fade(i,:))
            hold on;
            plot(1:size(mean_channels,2), mean(mean_channels(emg_dist.sess==i,:),1), 'Color', colors_fade(i,:), 'LineWidth', 2)
        end
        legend('sess 1', '', 'sess 2', '')
        legend boxoff
        title('AVG EMG activity','FontSize',my_font.title)
        xlim([0,length(mean_channels)+1])
        % ylim([0,0.8])
        xlabel('EMG channels','FontSize',my_font.xlabel)
        ylabel('AVG Normalized EMG')
        xticks(1:length(mean_channels))
        xticklabels(emg_locs_names)

    case 'chord_magnitude'
        plot_option = 1;
        measure = 'MD_efc';
        vararginoptions(varargin,{'plot_option','measure'});
        
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        subj = unique(data.sn);
        num_fingers = unique(data.num_fingers);
        
        % getting the values of measure:
        values = eval(['data.' measure]);
        
        % making the output of the model:
        out = [];
        for sn = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(sn)));

            for i = 1:length(sess)
                % calculating avg chord patterns:
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',['subj' sprintf('%.02d',subj(sn))],'plot_option',0,'normalize_channels',1,'sess',sess(i));
                n = get_num_active_fingers(chords);
                
                for j = 1:length(chords)
                    tmp.sn = subj(sn);
                    tmp.sess = sess(i);
                    tmp.num_fingers = n(j);
                    tmp.chordID = chords(j);
                    tmp.mag = vecnorm(pattern(chords==chords(j),:)')';
                    tmp.measure = values(data.chordID==chords(j) & data.sn==subj(sn) & data.sess==sess(i));
                    out = addstruct(out,tmp,'row',1);
                end
            end
        end
        
        % partial correlation between MD and magnitude considering
        % number of fingers effect:
        C  = [];
        mag_n = [];
        for sn = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(sn)));
            for i = 1:length(sess)
                X = out.mag(out.sess==sess(i) & out.sn==subj(sn));
                Y = out.measure(out.sess==sess(i) & out.sn==subj(sn));
                Z = make_design_matrix(out.chordID(out.sess==sess(i) & out.sn==subj(sn)),'n_fing');
                [rho,p,resX,resY] = partial_corr(X,Y,Z);
                
                tmp.sn = subj(sn);
                tmp.sess = sess(i);
                tmp.rho = rho;
                tmp.p = p;
                tmp.res_mag = {resX};
                mag_n = [mag_n ; resX];
                tmp.res_measure = {resY};
                C = addstruct(C,tmp,'row',1);
            end
        end
        
        % add mag_n and MD_n to out:
        out.mag_n = mag_n;

        if (plot_option)
            for sn = 1:length(subj)
                sess = unique(data.sess(data.sn==subj(sn)));
                for i = 1:length(sess)
                    figure;
                    hold on;
                    % single finger:
                    scatter_corr(out.mag(out.num_fingers==1 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==1 & out.sess==sess(i) & out.sn==subj(sn)), 'r', 'o'); hold on
                    % 3 finger:
                    scatter_corr(out.mag(out.num_fingers==3 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==3 & out.sess==sess(i) & out.sn==subj(sn)), 'b', 'filled')
                    % 5 finger:
                    scatter_corr(out.mag(out.num_fingers==5 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==5 & out.sess==sess(i) & out.sn==subj(sn)), 'k', 'filled')
                    title(sprintf('%s , sess %d',['subj' num2str(subj(sn))],sess(i)),'FontSize',my_font.title)
                    xlabel('Norm EMG','FontSize',my_font.xlabel)
                    ylabel('MD','FontSize',my_font.ylabel)
                    % xlim([0,3.6])
                    ylim([0,4])
                end
    
                for i = 1:length(sess)
                    figure;
                    hold on
                    scatter_corr(C.res_mag{C.sess==sess(i) & C.sn==subj(sn)}, C.res_measure{C.sess==sess(i) & C.sn==subj(sn)}, 'k', 'o')
                    title(sprintf('%s , sess %d',['subj' num2str(subj(sn))],sess(i)),'FontSize',my_font.title)
                    xlabel('Norm EMG (n regressed out)','FontSize',my_font.xlabel)
                    ylabel('MD (n regressed out)','FontSize',my_font.ylabel)
                    % xlim([-1.5,1.5])
                    ylim([-2,2])
                end
            end
        end
        varargout{1} = C;
        varargout{2} = out;
    
    case 'chord_coact'
        % coactivation of the muscles in chord EMGs:
        plot_option = 1;
        measure = 'MD_efc';
        vararginoptions(varargin,{'plot_option','measure'});
        
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        subj = unique(data.sn);
        num_fingers = unique(data.num_fingers);
        
        % getting the values of measure:
        values = eval(['data.' measure]);
        
        % making the output of the model:
        out = [];
        for sn = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(sn)));

            for i = 1:length(sess)
                % calculating avg chord patterns:
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',['subj' sprintf('%.02d',subj(sn))],'plot_option',0,'normalize_channels',1,'sess',sess(i));
                n = get_num_active_fingers(chords);
                
                for j = 1:length(chords)
                    tmp_pattern = pattern(chords==chords(j),:);
                    G = tmp_pattern' * tmp_pattern;
                    mag = trace(G);
                    coact = (sum(G(:)) - mag)/mag;
                    
                    tmp.sn = subj(sn);
                    tmp.sess = sess(i);
                    tmp.num_fingers = n(j);
                    tmp.chordID = chords(j);
                    tmp.mag = mag;
                    tmp.coact = coact;
                    tmp.measure = values(data.chordID==chords(j) & data.sn==subj(sn) & data.sess==sess(i));
                    out = addstruct(out,tmp,'row',1);
                end
            end
        end
        
        % partial correlation between MD and coact considering
        % number of fingers effect:
        C  = [];
        mag_n = [];
        for sn = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(sn)));
            for i = 1:length(sess)
                X = out.coact(out.sess==sess(i) & out.sn==subj(sn));
                Y = out.measure(out.sess==sess(i) & out.sn==subj(sn));
                Z = make_design_matrix(out.chordID(out.sess==sess(i) & out.sn==subj(sn)),'n_fing');
                [rho,p,resX,resY] = partial_corr(X,Y,Z);
                
                tmp.sn = subj(sn);
                tmp.sess = sess(i);
                tmp.rho = rho;
                tmp.p = p;
                tmp.res_mag = {resX};
                mag_n = [mag_n ; resX];
                tmp.res_measure = {resY};
                C = addstruct(C,tmp,'row',1);
            end
        end
        
        % add mag_n and MD_n to out:
        out.mag_n = mag_n;
        
        if (plot_option)
            for sn = 1:length(subj)
                sess = unique(data.sess(data.sn==subj(sn)));
                for i = 1:length(sess)
                    figure;
                    hold on;
                    % single finger:
                    scatter_corr(out.coact(out.num_fingers==1 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==1 & out.sess==sess(i) & out.sn==subj(sn)), 'r', 'o'); hold on
                    % 3 finger:
                    scatter_corr(out.coact(out.num_fingers==3 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==3 & out.sess==sess(i) & out.sn==subj(sn)), 'b', 'filled')
                    % 5 finger:
                    scatter_corr(out.coact(out.num_fingers==5 & out.sess==sess(i) & out.sn==subj(sn)), out.measure(out.num_fingers==5 & out.sess==sess(i) & out.sn==subj(sn)), 'k', 'filled')
                    title(sprintf('%s , sess %d',['subj' num2str(subj(sn))],sess(i)),'FontSize',my_font.title)
                    xlabel('Norm Coact EMG','FontSize',my_font.xlabel)
                    ylabel('MD','FontSize',my_font.ylabel)
                    % xlim([0,3.6])
                    ylim([0,4])
                end
    
                for i = 1:length(sess)
                    figure;
                    hold on
                    scatter_corr(C.res_mag{C.sess==sess(i) & C.sn==subj(sn)}, C.res_measure{C.sess==sess(i) & C.sn==subj(sn)}, 'k', 'o')
                    title(sprintf('%s , sess %d',['subj' num2str(subj(sn))],sess(i)),'FontSize',my_font.title)
                    xlabel('Norm Coact EMG (n regressed out)','FontSize',my_font.xlabel)
                    ylabel('MD (n regressed out)','FontSize',my_font.ylabel)
                    % xlim([-1.5,1.5])
                    ylim([-2,2])
                end
            end
        end
        varargout{1} = C;
        varargout{2} = out;

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
        d_type = 'project_to_nSphere';
        lambda = [];
        n_thresh = 5;
        sampling_option = 'whole_thresholded';
        plot_option = 1;
        vararginoptions(varargin,{'d_type','lambda','sampling_option','n_thresh','plot_option'})
        
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
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
            
            % gooz = natChord_analyze('EMG_prewhitening_matrix','plot_option',0);
            sess = unique(data.sess(data.sn == subjects(sn)));
            for j = 1:length(sess)
                % scaling factors:
                scales = get_emg_scales(subjects(sn),sess(j));

                % sigma = gooz.cov_res{gooz.sn==subjects(sn) & gooz.sess==sess(j)};
                
                emg_dist_sess = getrow(emg_dist,emg_dist.sess==sess(j));
                % normalizing the natural EMGs:
                for i = 1:length(emg_dist_sess.dist)
                    emg_dist_sess.dist{i} = (emg_dist_sess.dist{i} ./ scales); % * sigma^-(1/2): Ali did a prewhitening test
                end
    
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_names(sn,:),'plot_option',0,'normalize_channels',1,'sess',sess(j));
                n = get_num_active_fingers(chords);

                % pattern = gooz.pattern_prewhitened{gooz.sn==subjects(sn) & gooz.sess==sess(j)};
                
    
                % getting avg mean deviation of chords across sessions:
                chords_mean_dev = zeros(length(chords),1);
                for k = 1:length(chords)
                    row = data.sn==subjects(sn) & data.chordID==chords(k);
                    chords_mean_dev(k) = mean(data.MD_efc(row));
                end
                
                % container for the each session's dataframe:
                tmp = [];

                cnt = 1;
                for k = 1:length(chords)
                    d_avg = 0;
                    slope_avg = 0;
                    log_slope_avg = 0;
                    for i = 1:length(emg_dist_sess.partition(emg_dist_sess.sess==sess(j)))
                        % sorted distances from the natural dist:
                        d = get_d_from_natural(pattern(k,:)',emg_dist_sess.dist{emg_dist_sess.sess==sess(j) & emg_dist_sess.partition==i}, 'd_type', d_type, 'lambda',lambda);
                        
                        d_avg = d_avg + mean(d(1:n_thresh))/length(emg_dist_sess.partition(emg_dist_sess.sess==sess(j)));
                        slope_avg = slope_avg + linslope([d(1:n_thresh).^10,(1:n_thresh)'],'intercept',0)/length(emg_dist_sess.partition(emg_dist_sess.sess==sess(j))); %n_thresh/d(n_thresh)^10;
                        log_slope_avg = log_slope_avg + log(linslope([d(1:n_thresh).^10,(1:n_thresh)'],'intercept',0))/length(emg_dist_sess.partition(emg_dist_sess.sess==sess(j)));
                    end
                    % storing the information:
                    tmp.sn(cnt,1) = subjects(sn);
                    tmp.sess(cnt,1) = sess(j);
                    tmp.chordID(cnt,1) = chords(k);
                    tmp.num_fingers(cnt,1) = n(k);
                    tmp.MD(cnt,1) = chords_mean_dev(k);
                    tmp.thresh(cnt,1) = n_thresh;
                    tmp.d(cnt,1) = d_avg;
                    tmp.slope(cnt,1) = slope_avg; %n_thresh/d(n_thresh)^10;
                    tmp.log_slope(cnt,1) = log_slope_avg;

                    cnt = cnt+1;
                end
                C = addstruct(C,tmp,'row','force');
            end
        end

        % correlation of MD and log_slope:
        corr_struct = [];
        log_slope_n = [];
        MD_n = [];
        for sn = 1:length(subjects)
            for i = 1:length(unique(C.sess(C.sn==subjects(sn))))
                tmp = [];

                % correlation while regressing out the num finger effect:
                X = C.log_slope(C.sn==subjects(sn) & C.sess==i);
                Y = C.MD(C.sn==subjects(sn) & C.sess==i);
                Z = make_design_matrix(C.chordID(C.sn==subjects(sn) & C.sess==i), 'n_fing');
                [rho_n,p_n,res_log_slope_n,res_MD_n] = partial_corr(X,Y,Z);

                tmp.sn = subjects(sn);
                tmp.sess = i;
                tmp.rho_n = rho_n;
                tmp.p_n = p_n;
                log_slope_n = [log_slope_n;res_log_slope_n];
                MD_n = [MD_n;res_MD_n];
                corr_struct = addstruct(corr_struct,tmp,'row',1);
            end
        end

        % add log_slope_n and MD_n to C:
        C.log_slope_n = log_slope_n;
        C.MD_n = MD_n;

        if (plot_option)
            for sn = 1:length(subjects)
                for i = 1:length(unique(C.sess(C.sn==subjects(sn))))
                    figure;
                    % single finger:
                    hold on
                    scatter_corr(C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==1), C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==1), 'r', 'o')
                    % 3f chord:
                    scatter_corr(C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==3), C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==3), 'b', 'o')
                    % 5f chord:
                    scatter_corr(C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==5), C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==5), 'k', 'o')

                    title(sprintf('subj%d  , sess %d',subjects(sn),i),'FontSize',my_font.title)
                    xlabel('log(Slope (n/d))','FontSize',my_font.xlabel)
                    ylabel('MD','FontSize',my_font.ylabel)
                    % xlim([-5,25])
                    % ylim([0,4])
                end
                legend('1f','','','','','3f','','','','','5f','','','','')
                legend boxoff

                for i = 1:length(unique(C.sess(C.sn==subjects(sn))))
                    figure;
                    hold on
                    scatter_corr(C.log_slope_n(C.sn==subjects(sn) & C.sess==i), C.MD_n(C.sn==subjects(sn) & C.sess==i), 'k', 'o')
                    title(sprintf('subj%d  , sess %d',subjects(sn),i),'FontSize',my_font.title)
                    xlabel('log(Slope (n/d)) , n regressed out','FontSize',my_font.xlabel)
                    ylabel('MD , n regressed out','FontSize',my_font.ylabel)
                    xlim([-13,13])
                    ylim([-1.5,1.5])
                end
            end
        end

        varargout{1} = C;
        varargout{2} = corr_struct;


    case 'chord_distance_RDM'
        num_fingers = [];
        normalize_channels = 1;
        vararginoptions(varargin,{'subject_name','normalize_channels','num_fingers'})
        
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        data_all = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        sn_unique = unique(data.sn);

        C = [];

        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        % distance between muscle patterns:
        avg_pattern = 0;
        avg_d = 0;
        avg_mahalanobis = 0;
        for sn = 1:length(sn_unique)
            avg_pattern_sess = 0;
            avg_d_sess = 0;
            avg_mahalanobis_sess = 0;
            % temp struct:
            C_tmp = [];
            sess = unique(data.sess(data.sn==sn_unique(sn)));
            for i = 1:length(sess)
                % calculating avg chord patterns:
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',['subj' num2str(sn_unique(sn),'%.2d')],'plot_option',0,'normalize_channels',normalize_channels);
                
                % selected patterns of input chords:
                avg_pattern_sess = avg_pattern_sess + pattern/length(sess);

                % get Euclidean distance of patterns
                d = squareform(pdist(pattern));
                avg_d_sess = avg_d_sess + d/length(sess);

                % make the EMG pattern residuals required for Mahalanobis
                % distance:
                rows = data_all.sn==sn_unique(sn) & data_all.sess==i & data_all.trialCorr==1;
                emg_patterns_all = [data_all.emg_hold_avg_e1(rows), ...
                                    data_all.emg_hold_avg_e2(rows), ...
                                    data_all.emg_hold_avg_e3(rows), ...
                                    data_all.emg_hold_avg_e4(rows), ...
                                    data_all.emg_hold_avg_e5(rows), ...
                                    data_all.emg_hold_avg_f1(rows), ...
                                    data_all.emg_hold_avg_f2(rows), ...
                                    data_all.emg_hold_avg_f3(rows), ...
                                    data_all.emg_hold_avg_f4(rows), ...
                                    data_all.emg_hold_avg_f5(rows)];
                if (normalize_channels)
                    scales = get_emg_scales(sn_unique(sn),sess(i));
                    emg_patterns_all = emg_patterns_all ./ scales;
                end
                chords_all = data_all.chordID(rows);
                pattern_residuals = zeros(size(emg_patterns_all));
                for j = 1:length(chords)
                    idx_rows = chords_all==chords(j);
                    tmp_chord_patterns = emg_patterns_all(idx_rows,:);
                    pattern_residuals(idx_rows,:) = pattern(chords==chords(j),:) - tmp_chord_patterns;
                end
                
                % mahalanobis distance of patterns
                d_tmp = d_mahalanobis(pattern, cov(pattern_residuals));
                avg_mahalanobis_sess = avg_mahalanobis_sess + d_tmp/length(sess);
                % save the RDMs:
                C_tmp.sn = sn_unique(sn);
                C_tmp.sess = sess(i);
                C_tmp.d_mahalanobis{1} = d_tmp;
                C_tmp.chordID{1} = chords;
            end
            C = addstruct(C,C_tmp,'row','force');
            avg_pattern = avg_pattern + avg_pattern_sess/length(sn_unique);
            avg_d = avg_d + avg_d_sess/length(sn_unique);
            avg_mahalanobis = avg_mahalanobis + avg_mahalanobis_sess/length(sn_unique);
        end

        if ~isempty(num_fingers)
            n = get_num_active_fingers(chords);
            avg_pattern = avg_pattern(n==num_fingers,:);
            chords = chords(n==num_fingers);
            avg_d = avg_d(n==num_fingers,n==num_fingers);
            avg_mahalanobis = avg_mahalanobis(n==num_fingers,n==num_fingers);
            
            % bring flexions first:
            if (num_fingers==1)
                idx = [(6:10)';(1:5)'];
                chords = chords(idx);
                avg_pattern = avg_pattern(idx,:);

                tmpA = avg_d(6:10,6:10);
                tmpB = avg_d(1:5,1:5);
                tmpC = avg_d(1:5,6:10);
                tmpD = avg_d(6:10,1:5);
                avg_d(1:5,1:5) = tmpA;
                avg_d(1:5,6:10) = tmpD;
                avg_d(6:10,6:10) = tmpB;
                avg_d(6:10,1:5) = tmpC;

                tmpA = avg_mahalanobis(6:10,6:10);
                tmpB = avg_mahalanobis(1:5,1:5);
                tmpC = avg_mahalanobis(1:5,6:10);
                tmpD = avg_mahalanobis(6:10,1:5);
                avg_mahalanobis(1:5,1:5) = tmpA;
                avg_mahalanobis(1:5,6:10) = tmpD;
                avg_mahalanobis(6:10,6:10) = tmpB;
                avg_mahalanobis(6:10,1:5) = tmpC;
            end
        end

        figure;
        imagesc(avg_pattern); hold on;
        % clim([0, 1.6])
        colormap('parula')
        colorbar
        % plot settings:
        ax = gca;
        set(ax,'YTick',(1:size(chords,1)))
        set(ax,'YTickLabel',chords)
        set(ax,'XTick', (1:size(emg_locs_names,1)))
        set(ax,'XTickLabel',emg_locs_names)
        set(gca,'YDir','reverse')
        title(sprintf('Avg Muscle Patterns'))

        figure;
        imagesc(avg_d); hold on;
        % clim([0, 1.6])
        axis square
        colormap('magma')
        colorbar
        drawline(5.5,'dir','horz','linewidth',4)
        drawline(5.5,'dir','vert','linewidth',4)
        
        % plot settings:
        ax = gca;
        
        set(ax,'YTick',(1:size(chords,1)))
        set(ax,'YTickLabel',chords)
        
        set(ax,'XTick', (1:size(chords,1)))
        set(ax,'XTickLabel',chords)
        
        set(gca,'YDir','reverse')
        title(sprintf('AVG Euclidean Distance'))

        figure;
        imagesc(sqrt(abs(avg_mahalanobis))); hold on;
        % clim([0, 6.5])
        axis square
        colormap('magma')
        colorbar
        drawline(5.5,'dir','horz','linewidth',4)
        drawline(5.5,'dir','vert','linewidth',4)
        
        % plot settings:
        ax = gca;
        
        set(ax,'YTick',(1:size(chords,1)))
        set(ax,'YTickLabel',chords)
        
        set(ax,'XTick', (1:size(chords,1)))
        set(ax,'XTickLabel',chords)
        
        set(gca,'YDir','reverse')
        title(sprintf('AVG Mahalanobis Distance'))
        
        varargout{1} = C;
 
    case 'natural_distance_RDM'
        % handling input arguments:
        sampling_option = 'whole_sampled';
        vararginoptions(varargin,{'subject_name','sampling_option'});
        
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        sn_unique = unique(data.sn);

        % distance between muscle patterns:
        group_distance = 0;
        group_corr = 0;
        group_mahalanobis = 0;
        % loop on subjects:
        for sn = 1:length(sn_unique)
            % set natural EMG file name:
            file_name = fullfile(project_path, 'analysis', ['natChord_' ['subj' num2str(sn_unique(sn),'%.2d')] '_emg_natural_' sampling_option '.mat']);
            
            % loading natural EMG dists:
            emg_dist = load(file_name);
            emg_dist = emg_dist.emg_natural_dist;
            sess = unique(emg_dist.sess);
            partititons = unique(emg_dist.partition);
    
            % normalizing the natural EMGs:
            for i = 1:length(sess)
                % scaling factors:
                scales = get_emg_scales(sn_unique(sn),sess(i));
                for j = 1:length(partititons)
                    row = emg_dist.sess==sess(i) & emg_dist.partition==partititons(j);
                    emg_dist.dist{row} = emg_dist.dist{row} ./ scales;
                end
            end
            
            % Euclidean distance of EMG channels in natural:
            d_emg_sn = 0;
            corr_sn = 0;
            mahalanobis_sn = 0;
            for i = 1:length(sess)
                tmp = zeros(size(emg_dist.dist{1},2),size(emg_dist.dist{1},2));
                tmp_corr = 0;
                tmp_mahalanobis = 0;
                for j = 1:length(partititons)
                    row = emg_dist.sess==sess(i) & emg_dist.partition==partititons(j);
                    swapped_flex_extend = emg_dist.dist{row};
                    swapped_flex_extend = swapped_flex_extend(:,[6:10,1:5]);
                    % Euclidean distance of EMG Channels:
                    tmp = tmp + squareform(pdist(swapped_flex_extend'))/length(partititons);

                    % Correlation of EMG Channels:
                    tmp_corr = tmp_corr + cov(swapped_flex_extend)/length(partititons);

                    % tmp_mahalanobis = tmp_mahalanobis + d_mahalanobis(swapped_flex_extend')/length(partititons);;
                end
                d_emg_sn = d_emg_sn + tmp/length(sess);
                corr_sn = corr_sn + tmp_corr/length(sess);
                mahalanobis_sn = mahalanobis_sn + tmp_mahalanobis/length(sess);
            end

            group_distance = group_distance + d_emg_sn/length(sn_unique);
            group_corr = group_corr + corr_sn/length(sn_unique);
            group_mahalanobis = group_mahalanobis + mahalanobis_sn/length(sn_unique);
        end

        size(group_mahalanobis)
        
        % PLOT:
        figure;
        imagesc(1-group_corr); hold on;
        axis square
        colormap('magma')
        colorbar
        drawline(5.5,'dir','horz','linewidth',4)
        drawline(5.5,'dir','vert','linewidth',4)
        title('1 - rho of channels in Natural')
        % clim([0 1])
        
        % PLOT:
        figure;
        imagesc(group_distance); hold on;
        axis square
        colormap('magma')
        colorbar
        drawline(5.5,'dir','horz','linewidth',4)
        drawline(5.5,'dir','vert','linewidth',4)
        title('Distance Channels in Natural')

        % PLOT:
        figure;
        imagesc(group_mahalanobis); hold on;
        axis square
        colormap('magma')
        colorbar
        drawline(5.5,'dir','horz','linewidth',4)    
        drawline(5.5,'dir','vert','linewidth',4)
        title('Mahalanobis Distance Channels in Natural')


    case 'decoding'
        subject_name = 'subj01';
        vararginoptions(varargin,{'subject_name'})
        data = dload(fullfile(project_path,'analysis','natChord_all.tsv'));
        
        TN = kron((1:5)',ones(1,size(data.BN,1)/5));
        TN = TN(:);
        data = getrow(data,data.trialCorr==1);
        unique_chords = unique(data.chordID);
        TN(data.trialCorr ~= 1) = [];
        chords = num2str(data.chordID)-'0';
        
        % make design matirx:
        X = [data.emg_hold_avg_e1,data.emg_hold_avg_e2,data.emg_hold_avg_e3,data.emg_hold_avg_e4,data.emg_hold_avg_e5,...
             data.emg_hold_avg_f1,data.emg_hold_avg_f2,data.emg_hold_avg_f3,data.emg_hold_avg_f4,data.emg_hold_avg_f5];
        
        scales = [data.emg_baseline_e1,data.emg_baseline_e2,data.emg_baseline_e3,data.emg_baseline_e4,data.emg_baseline_e5,...
             data.emg_baseline_f1,data.emg_baseline_f2,data.emg_baseline_f3,data.emg_baseline_f4,data.emg_baseline_f5];
        
        X = X./scales;
        
        % make dependant vars:
        label = zeros(length(chords),1);
        Y = zeros(length(chords),10);
        for i = 1:size(chords,1)
            for j = 1:length(chords(i,:))
                if chords(i,j) == 1
                    Y(i,j) = 1;
                elseif chords(i,j) == 2
                    Y(i,j+5) = 1;
                end
            end
            label(i) = find(unique_chords==data.chordID(i));
        end

        % train and test single finger:
        X_train = X(data.num_fingers==1 & TN~=5,:);
        Y_train = Y(data.num_fingers==1 & TN~=5,:);
        X_test = X(data.num_fingers==1 & TN==5,:);
        Y_test = Y(data.num_fingers==1 & TN==5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);
        
        % train single finger, test multi finger:
        X_train = X(data.num_fingers==1,:);
        Y_train = Y(data.num_fingers==1,:);
        X_test = X(data.num_fingers~=1,:);
        Y_test = Y(data.num_fingers~=1,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);

        % train and test overall:
        X_train = X(TN ~= 5,:);
        Y_train = Y(TN ~= 5,:);
        X_test = X(TN == 5,:);
        Y_test = Y(TN == 5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);

        % train and test overall with labels:
        X_train = X(TN ~= 5,:);
        Y_train = label(TN ~= 5,:);
        X_test = X(TN == 5,:);
        Y_test = label(TN == 5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test,'label',1);
        fprintf('label: overall acc = %.4f\n',C.acc_test);

        % train and test single finger with labels:
        X_train = X(data.num_fingers==1 & TN~=5,:);
        Y_train = label(data.num_fingers==1 & TN~=5,:);
        X_test = X(data.num_fingers==1 & TN==5,:);
        Y_test = label(data.num_fingers==1 & TN==5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test,'label',1);
        fprintf('label: 1f acc = %.4f\n',C.acc_test);

        % predict EMG from chordID - single finger:
        X_train = Y(data.num_fingers==1 & TN ~= 5,:);
        Y_train = X(data.num_fingers==1 & TN ~= 5,:);
        X_test = Y(data.num_fingers==1 & TN == 5,:);
        Y_test = X(data.num_fingers==1 & TN == 5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);
        C.r

        % predict EMG from chordID - chord from 1f:
        X_train = Y(data.num_fingers==1,:);
        Y_train = X(data.num_fingers==1,:);
        X_test = Y(data.num_fingers~=1,:);
        Y_test = X(data.num_fingers~=1,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);
        C.r

        % predict EMG from chordID - overall:
        X_train = Y(TN ~= 5,:);
        Y_train = X(TN ~= 5,:);
        X_test = Y(TN == 5,:);
        Y_test = X(TN == 5,:);
        C = linear_classify(X_train,Y_train,X_test,Y_test);
        C.r

        varargout{1} = C;

    case 'test_nSphere_magnitude'
        C = natChord_analyze('nSphere_model','d_type','project_to_nSphere','n_thresh',1,'sampling_option','whole_thresholded','plot_option',0);
        subjects = unique(C.sn);

        out = [];
        for sn = 1:length(subjects)
            mag = natChord_analyze('chord_magnitude_difficulty_model','subject_name',['subj' sprintf('%02d',subjects(sn))],'plot_option',0);
            for i = 1:length(unique(C.sess))
                % 1f:
                MD_1f = C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==1);
                log_1f = C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==1);
                mag_1f = mag.mag(mag.sess==i & mag.num_fingers==1);

                % 3f:
                MD_3f = C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==3);
                log_3f = C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==3);
                mag_3f = mag.mag(mag.sess==i & mag.num_fingers==3);

                % 5f:
                MD_5f = C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers==5);
                log_5f = C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers==5);
                mag_5f = mag.mag(mag.sess==i & mag.num_fingers==5);

                % all chords:
                MD_chords = C.MD(C.sn==subjects(sn) & C.sess==i & C.num_fingers>1);
                log_chords = C.log_slope(C.sn==subjects(sn) & C.sess==i & C.num_fingers>1);
                mag_chords = mag.mag(mag.sess==i & mag.num_fingers>1);

                % merged model:
                mdl_3f = fitlm([log_3f,mag_3f], MD_3f);
                mdl_5f = fitlm([log_5f,mag_5f], MD_5f);
                mdl_chords = fitlm([log_chords,mag_chords], MD_chords);
                
                tmp.sn = subjects(sn);
                tmp.sess = i;

                tmp.mag_corr_3f = corr(mag_3f,MD_3f);
                tmp.mag_corr_5f = corr(mag_5f,MD_5f);
                tmp.mag_corr_chords = corr(mag_chords,MD_chords);

                tmp.likelihood_corr_3f = corr(log_3f,MD_3f);
                tmp.likelihood_corr_5f = corr(log_5f,MD_5f);
                tmp.likelihood_corr_chords = corr(log_chords,MD_chords);

                tmp.merged_corr_3f = corr(mdl_3f.Fitted,MD_3f);
                tmp.merged_corr_5f = corr(mdl_5f.Fitted,MD_5f);
                tmp.merged_corr_chords = corr(mdl_chords.Fitted,MD_chords);

                out = addstruct(out,tmp,'row',1);
            end
        end

        varargout{1} = out;

    case 'component_angle_model'
        % handling input arguments:
        sampling_option = 'whole_sampled';
        n_components = 5;
        plot_option = 1;
        vararginoptions(varargin,{'sampling_option','plot_option','n_components'})
        
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        sess = unique(data.sess);
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

            for j = 1:length(sess)
                % scaling factors:
                scales = get_emg_scales(subjects(sn),sess(j));
                
                % normalizing the natural EMGs:
                for i = 1:length(emg_dist.dist)
                    emg_dist.dist{i} = emg_dist.dist{i} ./ scales;
                end

                for i = 1:length(emg_dist.dist)
                    tmp = transform_dist(emg_dist.dist{i},n_components);
                    emg_dist.dist{i} = tmp.transformed_dist{1};
                    emg_dist.transform_mat(i,1) = tmp.transform_mat;
                    emg_dist.mdl(i,1) = tmp.mdl;
                end
            
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',subject_names(sn,:),'plot_option',0,'normalize_channels',1,'sess',sess(j));
                n = get_num_active_fingers(chords);
                
                % getting avg mean deviation of chords across sessions:
                chords_mean_dev = zeros(length(chords),1);
                for k = 1:length(chords)
                    row = data.sn==subjects(sn) & data.chordID==chords(k);
                    chords_mean_dev(k) = data.MD_efc(row & data.sess==sess(j));
                end
                
                for k = 1:length(chords)
                    tmp = [];
                    gamma_avg = 0;
                    norm_avg = 0;
                    for i = 1:length(emg_dist.partition(emg_dist.sess==sess(j)))
                        ica_components = emg_dist.transform_mat{i};

                        % find the angle with all components:
                        gamma_tmp = zeros(size(ica_components,2),1);
                        norm_tmp = zeros(size(ica_components,2),1);
                        for i_component = 1:size(ica_components,2)
                            gamma_tmp(i_component) = abs(cos_angle(pattern(k,:)',ica_components(:,i_component)));
                            
                            vec1 = pattern(k,:)'/norm(pattern(k,:)');
                            vec2 = ica_components(:,i_component);
                            norm_tmp(i_component) = norm(vec1-vec2);
                        end
                        gamma_avg = gamma_avg + max(gamma_tmp)/length(emg_dist.partition(emg_dist.sess==sess(j)));
                        norm_avg = norm_avg + min(norm_tmp)/length(emg_dist.partition(emg_dist.sess==sess(j)));
                    end

                    % storing the values:
                    tmp.sn = subjects(sn);
                    tmp.sess = sess(j);
                    tmp.chordID = chords(k);
                    tmp.num_fingers = n(k);
                    tmp.MD = chords_mean_dev(k);
                    tmp.gamma = gamma_avg;
                    tmp.norm = norm_avg;

                    C = addstruct(C,tmp,'row','force');
                end
            end
        end

        corr_struct = [];
        for sn = 1:length(subjects)
            for j = 1:length(sess)
                tmp = [];
                tmp.sn = subjects(sn);
                tmp.sess = j;
                tmp.rho_1f = corr(C.gamma(C.num_fingers==1 & C.sn==subjects(sn) & C.sess==j),C.MD(C.num_fingers==1 & C.sn==subjects(sn) & C.sess==j));
                tmp.rho_3f = corr(C.gamma(C.num_fingers==3 & C.sn==subjects(sn) & C.sess==j),C.MD(C.num_fingers==3 & C.sn==subjects(sn) & C.sess==j));
                tmp.rho_5f = corr(C.gamma(C.num_fingers==5 & C.sn==subjects(sn) & C.sess==j),C.MD(C.num_fingers==5 & C.sn==subjects(sn) & C.sess==j));
            
                corr_struct = addstruct(corr_struct,tmp,'row',1);
            end
        end

        if plot_option
            for sn = 1:length(subjects)
                for j = 1:length(sess)
                    X = C.gamma(C.sn==subjects(sn) & C.sess==j);
                    Y = C.MD(C.sn==subjects(sn) & C.sess==j);
                    Z = make_design_matrix(C.chordID(C.sn==subjects(sn) & C.sess==j),'n_fing');
                    [rho,p,resX,resY] = partial_corr(X,Y,Z);
                    figure;
                    hold on;
                    scatter_corr(resX,resY,'k','o')
                    xlabel('max cosine angle, n regressed out')
                    ylabel('MD, n regressed out')
                end
            end
        end

        varargout{1} = corr_struct;
        varargout{2} = C;

    case 'model_testing_all'
        % handling input arguments:
        measure = 'MD_efc';
        model_names = {'n_fing','n_fing+nSphere','n_fing+magnitude','n_fing+magnitude+nSphere'};
        vararginoptions(varargin,{'chords','measure','model_names'})
        
        % loading data:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(1:sum(data.sn==1 & data.sess==1));
        subj = unique(data.sn);
        sess = unique(data.sess);
        
        % getting the values of measure:
        values_tmp = eval(['data.' measure]);
        
        % getting the average of sessions for every subj:
        values = zeros(length(chords),length(subj));
        for i = 1:length(subj)
            % avg with considering nan values since subjects might have
            % missed all 5 repetitions in one session:
            tmp = [];
            for j = 1:length(sess)
                tmp = [tmp, values_tmp(data.sess==sess(j) & data.sn==subj(i))];
            end
            values(:,i) = mean(tmp,2,'omitmissing');
        end
        
        % modelling the difficulty for all chords.
        C = [];
        % loop on subjects and model testing witihin subjects:
        for sn = 1:length(subj)
            % values of subject, Nx1 vector:
            y = values(:,sn);
            
            % loop on models to be tested:
            for i_mdl = 1:length(model_names)
                % getting design matrix for model:
                X = make_design_matrix(chords,model_names{i_mdl},'sn',subj(sn));

                % check design matrix's Rank:
                is_full_rank = rank(X) == size(X,2);

                % training the model:
                [B,STATS] = svd_linregress(y,X);

                % testing the model:
                y_pred = X*B;
                r = corr(y_pred,y);
                SSR = sum((y_pred-y).^2);
                SST = sum((y-mean(y)).^2);
                r2 = 1 - SSR/SST;

                % storing the results:
                C_tmp.model = model_names(i_mdl);
                C_tmp.is_full_rank = is_full_rank;
                C_tmp.B = {B};
                C_tmp.stats = {STATS};
                C_tmp.r = r;
                C_tmp.r2 = r2;

                C = addstruct(C,C_tmp,'row',1);
            end
        end

        % PLOT - regression results:
        figure;
        ax1 = axes('Units', 'centimeters', 'Position', [2 2 3.5 3],'Box','off');
        for j = 1:length(model_names)
            % getting cross validated r:
            r = C.r(strcmp(C.model,model_names{j}));
            
            r_avg(j) = mean(r);
            r_sem(j) = std(r)/sqrt(length(r));
        end
        drawline(1,'dir','horz','color',[0.7 0.7 0.7],'lim',[0,length(model_names)+1],'linestyle',':'); hold on;
        plot(1:length(model_names),r_avg,'LineWidth',2,'Color',[0.1 0.1 0.1,0.1]);
        errorbar(1:length(model_names),r_avg,r_sem,'LineStyle','none','Color',[0.1 0.1 0.1],'CapSize',0)
        scatter(1:length(model_names),r_avg,15,'filled','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1]);
        box off
        h = gca;
        h.YTick = 0:0.25:1;
        h.XTick = 1:length(model_names);
        h.XTickLabel = cellfun(@(x) replace(x,'_',' '),model_names,'uniformoutput',false);
        h.XAxis.FontSize = my_font.tick_label;
        h.YAxis.FontSize = my_font.tick_label;
        ylabel('r','FontSize',my_font.ylabel)
        ylim([0, 1.05])
        xlim([0.5,length(model_names)+0.5])
        fontname("Arial")
        
        varargout{1} = C;

    case 'repetition_effect'
        measure = 'MD';
        vararginoptions(varargin,{'chords','measure'})
        
        data = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        chords = unique(data.chordID);
        n = get_num_active_fingers(chords);

        % getting the values of measure:
        values = eval(['data.' measure]);
        values(values==-1) = NaN;
        
        % putting trials in rows:
        n_fing = reshape(data.num_fingers,5,[]);
        sess = reshape(data.sess,5,[]);
        values = reshape(values,5,[]);
        subj = reshape(data.sn,5,[]);
        chordID = reshape(data.chordID,5,[]);
        repetitions = 5;

        % getting averages within each session and n_fing:
        n_fing = n_fing(1,:);
        n_fing_unique = unique(n_fing);
        chordID = chordID(1,:);
        sess = sess(1,:);
        subj = subj(1,:);
        subj_unique = unique(subj);
        C = [];
        % loop on num fingers:
        for i = 1:length(unique(n_fing_unique))
            for sn = 1:length(subj_unique)
                subj_sess = unique(sess(subj==subj_unique(sn)));
                for j = 1:length(subj_sess)
                    % variable to hold the values of all the repetitions of
                    % the chords. Chords are stacked on top of each other.
                    % Meaning that each row of the tmp_container will be
                    % a specific chord and columns are all the repetitions.
                    % Therefore, number of columns will be 'repetition*number of times chord was visited in the session'
                    tmp_container = [];
                    chords_tmp = chords(n==n_fing_unique(i));
                    for k = 1:length(chords_tmp)
                        % getting all the data that chord was repeated. the
                        % columns are the times that the chord was
                        % repeated in each session:
                        values_tmp = values(:, subj==subj_unique(sn) & sess==j & n_fing==n_fing_unique(i) & chordID==chords_tmp(k));
                        num_visited = size(values_tmp,2);
    
                        % this condition is because two of the chords were
                        % presented four times in each session (by mistake):
                        if (num_visited==4)
                            values_tmp = (values_tmp(:,1:2)+values_tmp(:,3:4))*0.5;
                        end
                        tmp_container = [tmp_container ; values_tmp(:)'];
                    end
                    
                    num_visited = size(tmp_container,2)/repetitions;
                    for k = 1:num_visited
                        tmp = [];
                        tmp.sn = subj_unique(sn);
                        tmp.sess = j;
                        tmp.num_fingers = n_fing_unique(i);
                        % the time that the chord was visited during the
                        % sessions.
                        tmp.visit = k;
                        tmp.value_subj(1,:) = mean(tmp_container(:,(k-1)*repetitions+1:k*repetitions),1,'omitmissing'); 
                        C = addstruct(C,tmp,'row',1);
                    end
                end
            end
        end
        
        % PLOT - repetition trends across visits and sessions:
        sess = unique(C.sess);
        for j = 1:length(sess)
            figure;
            ax1 = axes('Units', 'centimeters', 'Position', [2 2 2.4 5],'Box','off');
            visits = unique(C.visit);
            offset_size = 5;
            x_offset = 0:offset_size:5*(length(visits)-1);
            for i = 1:length(n_fing_unique)
                for k = 1:length(visits)
                    row = C.num_fingers==n_fing_unique(i) & C.visit==k & C.sess==j;
                    plot((1:repetitions)+x_offset(k), mean(C.value_subj(row, :),1),'Color',colors_blue(n_fing_unique(i),:),'LineWidth',1); hold on;
                    % errorbar((1:repetitions)+x_offset(k)+(j-1)*offset_size*2, mean(C.value_subj(row, :),1), mean(C.sem(C.num_fingers==i & C.sess==j, :),1), 'CapSize', 0, 'Color', colors_blue(n_fing_unique(i),:));
                    scatter((1:repetitions)+x_offset(k), mean(C.value_subj(row, :),1), 10,'MarkerFaceColor',colors_blue(n_fing_unique(i),:),'MarkerEdgeColor',colors_blue(n_fing_unique(i),:))
                end
            end
            box off
            h = gca;
            h.YTick = 100:150:600; % RT
            % h.YTick = 0:1000:3000; % MT
            % h.YTick = 0.5:1:2.5; % MD
            h.XTick = 5*(1:length(visits)) - 2;
            xlabel('visit','FontSize',my_font.xlabel)
            h.XTickLabel = {'1','2'};
            h.XAxis.FontSize = my_font.tick_label;
            h.YAxis.FontSize = my_font.tick_label;
            ylabel(measure,'FontSize',my_font.ylabel)
            % ylabel([measure ' [ms]'],'FontSize',my_font.ylabel)
            % ylim([0.3, 3]) % MD
            ylim([0, 650]) % RT
            % ylim([0, 3200]) % MT
            xlim([0,11])
            % title('Repetition Effect','FontSize',my_font.title)
            fontname("Arial")
        end

        varargout{1} = C;

    
    case 'likelihood_reliability'
        measure = 'log_slope';
        centered = 1;
        vararginoptions(varargin,{'measure','centered'})

        [~,C] = natChord_analyze('nSphere_model','plot_option',0);

        % getting the values of measure:
        values = eval(['C.' measure]);
        
        % reliability estimation:
        [v_g, v_gs, v_gse] = reliability_var(values, C.sn, C.sess, ...
            'cond_vec', C.num_fingers, 'centered', centered);
        
        % plot:
        y = [];
        figure;
        ax1 = axes('Units','centimeters', 'Position', [2 2 4 4],'Box','off');
        for i = 1:length(unique(C.num_fingers))
            y(i,:) = [v_g{i}/v_gse{i} (v_gs{i}-v_g{i})/v_gse{i} (v_gse{i}-v_gs{i})/v_gse{i}];
            b = bar(i,y(i,:),'stacked','FaceColor','flat','BarWidth',0.8);
            b(1).CData = colors_blue(5,:);   % global var
            b(2).CData = colors_red(3,:);  % subj var
            b(3).CData = [0.8 0.8 0.8];  % noise var
            b(1).EdgeColor = [1 1 1];
            b(2).EdgeColor = [1 1 1];
            b(3).EdgeColor = [1 1 1];
            b(1).LineWidth = 1;
            b(2).LineWidth = 1;
            b(3).LineWidth = 1;
            hold on
        end
        drawline(1,'dir','horz','color',[0.7 0.7 0.7])
        % drawline(0,'dir','horz','color',[0.7 0.7 0.7])
        h = gca;
        h.XTick = 1:length(unique(C.num_fingers));
        h.XTickLabel = num2str(unique(C.num_fingers));
        h.YTick = 0:0.2:1;
        box off;
        % title([measure ' Reliability'],'FontSize',my_font.title)
        xlabel('num fingers','FontSize',my_font.xlabel)
        ylabel('percent variance','FontSize',my_font.ylabel)
        % lgd = legend('global','subject','noise');
        % legend boxoff
        % fontsize(lgd,6,'points')
        ylim([0,1.1])
        xlim([0.3,3.7])
        h.XAxis.FontSize = my_font.tick_label;
        h.YAxis.FontSize = my_font.tick_label;
        fontname("Arial")

        varargout{1} = [v_g, v_gs, v_gse];

    case 'model_testing_all_efc1'
        % handling input arguments:
        measure = 'MD';
        sess = [3,4];
        model_names = {'n_fing','emg_additive_avg','n_fing+transition','n_fing+additive','n_fing+all_2fing','n_fing+additive+all_2fing',... % behavioural models
                       'n_fing+force_avg','n_fing+force_2fing','n_fing+force_avg+force_2fing',... % force models
                       'n_fing+emg_additive_avg','n_fing+emg_2channel_avg','n_fing+emg_additive_avg+emg_2channel_avg',... % EMG models
                       'n_fing+nSphere_avg','n_fing+magnitude_avg','n_fing+nSphere_avg+magnitude_avg'}; % inference models
        % model_names = {'n_fing','n_fing+transition','n_fing+additive','n_fing+symmetries','n_fing+all_2fing','n_fing+additive+all_2fing',... % behavioural models
        %                'n_fing+force_avg','n_fing+force_2fing','n_fing+force_avg+force_2fing',... % force models
        %                'n_fing+emg_additive_avg','n_fing+emg_2channel_avg','n_fing+emg_additive_avg+emg_2channel_avg',... % EMG models
        %                'n_fing+nSphere_avg','n_fing+magnitude_avg','n_fing+coact_avg','n_fing+magnitude_avg+coact_avg','n_fing+nSphere_avg+magnitude_avg','n_fing+nSphere_avg+coact_avg','n_fing+nSphere_avg+coact_avg+magnitude_avg'}; % inference models
        base_models = {'n_fing','transition','additive','all_2fing','force_avg','force_2fing','emg_additive_avg','emg_2channel_avg',...
                       'nSphere_avg','magnitude_avg'};
        % base_models = {'n_fing','transition','additive','symmetries','all_2fing','force_avg','force_2fing','emg_additive_avg','emg_2channel_avg',...
        %                'nSphere_avg','magnitude_avg', 'coact_avg'};
        vararginoptions(varargin,{'chords','measure','model_names'})
        
        % loading data:
        chords_natChord = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords_natChord = chords_natChord.chordID(chords_natChord.sn==1 & chords_natChord.sess==1);
        data = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        data = getrow(data,ismember(data.chordID,chords_natChord));
        chords = data.chordID(data.sn==1 & data.sess==1);
        subj = unique(data.sn);
        
        % getting the values of measure:
        values_tmp = eval(['data.' measure]);
        
        % getting the average of sessions for every subj:
        values = zeros(length(chords),length(subj));
        for i = 1:length(subj)
            % avg with considering nan values since subjects might have
            % missed all 5 repetitions in one session:
            values(:,i) = mean([values_tmp(data.sess==sess(1) & data.sn==subj(i)),values_tmp(data.sess==sess(2) & data.sn==subj(i))],2,'omitmissing');
        end

        % modelling the difficulty for chords
        C = [];
        % loop on subjects and regression with leave-one-out:
        for sn = 1:length(subj)
            % values of 'in' subjects, Nx1 vector:
            y_train = values(:,setdiff(1:length(subj),sn));
            y_train = mean(y_train,2);

            % avg of 'out' subject:
            y_test = values(:,sn);

            % loop on models to be tested:
            for i_mdl = 1:length(model_names)
                % getting design matrix for model:
                X = make_design_matrix(chords,model_names{i_mdl});
                % X = [ones(size(X,1),1) X];

                % check design matrix's Rank:
                is_full_rank = rank(X) == size(X,2);

                % training the model:
                % [B,STATS] = linregress(y_train,X,'intercept',0);
                [B,STATS] = svd_linregress(y_train,X);

                % testing the model:
                X_test = make_design_matrix(chords,model_names{i_mdl});
                % X_test = [ones(size(X_test,1),1) X_test];
                y_pred = X_test * B;
                r = corr(y_pred,y_test);
                SSR = sum((y_pred-y_test).^2);
                SST = sum((y_test-mean(y_test)).^2);
                r2 = 1 - SSR/SST;

                % storing the results:
                C_tmp.sn_out = sn;
                C_tmp.model = model_names(i_mdl);
                C_tmp.is_full_rank = is_full_rank;
                C_tmp.B = {B};
                C_tmp.stats = {STATS};
                C_tmp.r = r;
                C_tmp.r2 = r2;

                C = addstruct(C,C_tmp,'row',1);
            end
        end

        [R,R2] = crossval_reliability(values);
        noise_ceil = mean(R);

        C.R_ceil = repmat(R,length(model_names),1);
        C.R2_ceil = repmat(R2,length(model_names),1);

        % analysis dataframe:
        df = [];
        df.sn_out = C.sn_out;
        df.r = C.r;
        df.r2 = C.r2;
        df.r_ceil = C.R_ceil;
        df.r2_ceil = C.R2_ceil;
        % make flat strcture for models:
        for i = 1:length(base_models)
            df.(base_models{i}) = zeros(length(df.sn_out),1);
        end
        % fill the model name rows in the df:
        for i = 1:length(model_names)
            rows = strcmp(C.model,model_names{i});
            split_model_names = strsplit(model_names{i},'+');
            for j = 1:length(split_model_names)
                df.(split_model_names{j})(rows) = 1;
            end
        end
        dsave(fullfile(project_path,'analysis',['emg_models_' measure '.tsv']),df);
        save(fullfile(project_path,'analysis',['emg_models_details_' measure '.mat']),'C')
        
        % stats between models:
        stats = [];
        for i = 1:length(model_names)-1
            r1 = C.r(strcmp(C.model,model_names{i}));
            for j = i+1:length(model_names)
                r2 = C.r(strcmp(C.model,model_names{j}));
                % paired t-test, one-tail r2>r1:
                [t,p] = ttest(r2,r1,1,'paired');
                tmp.models = {model_names{i},model_names{j}};
                tmp.t = t;
                tmp.p = p;
                stats = addstruct(stats,tmp,'row',1);
            end
        end
        varargout{1} = C;
        varargout{2} = stats;


    case 'model_testing_learning_slope_efc1'
        % handling input arguments:
        measure = 'MD';
        model_names = {'n_fing','n_fing+magnitude_avg','n_fing+magnitude_avg+nSphere_avg','n_fing+transition','n_fing+additive+2fing_adj'};
        vararginoptions(varargin,{'chords','measure','model_names'})
        
        % loading data:
        chords_natChord = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords_natChord = chords_natChord.chordID(chords_natChord.sn==2);
        data = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        data = getrow(data,ismember(data.chordID,chords_natChord));
        chords = data.chordID(data.sn==1 & data.sess==1);
        subj = unique(data.sn);
        
        % getting the values of measure:
        values_tmp = eval(['data.' measure]);
        
        % getting the average of sessions for every subj:
        values = zeros(length(chords),length(subj));
        for i = 1:length(subj)
            tmp = [values_tmp(data.sess==1 & data.sn==subj(i)), values_tmp(data.sess==2 & data.sn==subj(i)), ...
                   values_tmp(data.sess==3 & data.sn==subj(i)), values_tmp(data.sess==4 & data.sn==subj(i))];

            for j = 1:size(tmp,1)
                tmp_x = (1:4)';
                tmp_y = tmp(j,:)';
                nan_idx = isnan(tmp_y);
                tmp_y(nan_idx) = [];
                tmp_x(nan_idx) = [];
                tmp_slope = linslope([tmp_x , tmp_y]);

                values(j,i) = tmp_slope;
            end
        end

        % modelling the learning slope for all chords.
        C = [];
        % loop on subjects and regression with leave-one-out:
        for sn = 1:length(subj)
            % values of 'in' subjects, Nx1 vector:
            y_train = values(:,setdiff(1:length(subj),sn));
            y_train = mean(y_train,2);

            % avg of 'out' subject:
            y_test = values(:,sn);

            % loop on models to be tested:
            for i_mdl = 1:length(model_names)
                % getting design matrix for model:
                X = make_design_matrix(chords,model_names{i_mdl});

                % check design matrix's Rank:
                is_full_rank = rank(X) == size(X,2);

                % training the model:
                % [B,STATS] = linregress(y_train,X,'intercept',0);
                [B,STATS] = svd_linregress(y_train,X);

                % testing the model:
                X_test = make_design_matrix(chords,model_names{i_mdl});
                y_pred = X_test * B;
                r = corr(y_pred,y_test);
                SSR = sum((y_pred-y_test).^2);
                SST = sum((y_test-mean(y_test)).^2);
                r2 = 1 - SSR/SST;

                % storing the results:
                C_tmp.sn_out = sn;
                C_tmp.model = model_names(i_mdl);
                C_tmp.is_full_rank = is_full_rank;
                C_tmp.B = {B};
                C_tmp.stats = {STATS};
                C_tmp.r = r;
                C_tmp.r2 = r2;

                C = addstruct(C,C_tmp,'row',1);
            end
        end
        
        % stats between models:
        stats = [];
        for i = 1:length(model_names)-1
            r1 = C.r(strcmp(C.model,model_names{i}));
            for j = i+1:length(model_names)
                tmp = [];
                r2 = C.r(strcmp(C.model,model_names{j}));
                % paired t-test, one-tail r2>r1:
                [t,p] = ttest(r2,r1,1,'paired');
                tmp.models = {model_names{i},model_names{j}};
                tmp.t = t;
                tmp.p = p;
                stats = addstruct(stats,tmp,'row',1);
            end
        end


        % getting noise ceiling:
        % [~,corr_struct] = efc1_analyze('selected_chords_reliability','blocks',[(sess(1)-1)*12+1 sess(2)*12],'chords',chords,'plot_option',0);
        % if (strcmp(measure,'MD'))
        %     noise_ceil = mean(corr_struct.MD);
        % elseif (strcmp(measure,'MT'))
        %     noise_ceil = mean(corr_struct.MT);
        % else
        %     noise_ceil = mean(corr_struct.RT);
        % end
        noise_ceil = 0.8685;

        for i = 1:length(model_names)
            r = C.r(strcmp(C.model,model_names{i}));
            fprintf('ttest: model %s different from noise ceiling:\n',model_names{i})
            ttest(r-noise_ceil,[],2,'onesample')
        end

        % PLOT - regression results:
        figure;
        ax1 = axes('Units', 'centimeters', 'Position', [2 2 3.5 3],'Box','off');
        for j = 1:length(model_names)
            % getting cross validated r:
            r = C.r(strcmp(C.model,model_names{j}));
            
            r_avg(j) = mean(r);
            r_sem(j) = std(r)/sqrt(length(r));
        end
        drawline(noise_ceil,'dir','horz','color',[0.7 0.7 0.7],'lim',[0,length(model_names)+1],'linestyle',':'); hold on;
        plot(1:length(model_names),r_avg,'LineWidth',2,'Color',[0.1 0.1 0.1,0.1]);
        errorbar(1:length(model_names),r_avg,r_sem,'LineStyle','none','Color',[0.1 0.1 0.1],'CapSize',0)
        scatter(1:length(model_names),r_avg,15,'filled','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1]);
        box off
        h = gca;
        h.YTick = 0:0.25:1;
        h.XTick = 1:length(model_names);
        h.XTickLabel = cellfun(@(x) replace(x,'_',' '),model_names,'uniformoutput',false);
        h.XAxis.FontSize = my_font.tick_label;
        h.YAxis.FontSize = my_font.tick_label;
        ylabel('r','FontSize',my_font.ylabel)
        % ylim([0, 0.2])
        xlim([0.5,length(model_names)+0.5])
        fontname("Arial")

        varargout{1} = C;
        varargout{2} = stats;

    case 'chord_natural_EMG_patterns'        
        % load data:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);

        emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];

        % get group avg chord patterns:
        patterns = make_design_matrix(chords,'emg_additive_avg','sess',1);
        rows = [1,6,10,11,23,25,26:4:68];
        patterns = patterns(rows,:);
        
        % loading natural EMG dists:
        subject_name = 'subj08';
        emg_dist = load(fullfile(project_path,'analysis',['natChord_' subject_name '_emg_natural_whole_sampled.mat']));
        emg_dist = emg_dist.emg_natural_dist;
        emg_dist = getrow(emg_dist,emg_dist.sess==1);

        % scaling factors:
        scales = get_emg_scales(str2double(subject_name(end-1:end)),1);

        % normalizing the natural EMGs:
        for i = 1:length(emg_dist.dist)
            emg_dist.dist{i} = emg_dist.dist{i} ./ scales;
        end

        
        % PLOT - Chords:
        figure;
        % plotting the chord EMG patterns:
        h = pcolor([[patterns, zeros(size(patterns,1),1)] ; zeros(1,size(patterns,2)+1)]);
        colorbar
        % plot settings:
        set(h,'EdgeColor','none')
        % plot settings:
        ax = gca;
        set(ax,'YTick',(1:size(patterns,1))+0.5)
        set(ax,'YTickLabel',chords)
        set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
        set(ax,'XTickLabel',emg_locs_names)
        set(gca,'YDir','reverse')
        
        
        % PLOT - Natural:
        i = 1;
        n = 100;
        matrix = emg_dist.dist{i}(2001:4800,:);
        % Reshape the matrix
        [~, numCols] = size(matrix);
        
        reshaped_matrix = reshape(matrix, n, [], numCols);
        
        % Take the mean along the first dimension
        averaged_matrix = mean(reshaped_matrix, 1);
        averaged_matrix = squeeze(averaged_matrix)'; % Remove singleton dimensions
        
        % swap channel locs in matrix:
        swap = averaged_matrix(6:10,:);
        averaged_matrix(6:10,:) = averaged_matrix(1:5,:);
        averaged_matrix(1:5,:) = swap;

        % plotting the natural EMG patterns:
        figure;
        h  = pcolor([[averaged_matrix, zeros(size(averaged_matrix,1),1)] ; zeros(1,size(averaged_matrix,2)+1)]);
        colorbar
        % plot settings:
        set(h, 'EdgeColor', [0.2,0.2,0.2]);
        % ax = gca;
        % set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
        % set(ax,'XTickLabel',emg_locs_names)
        set(gca,'YDir','reverse')
        title('natural')
        % title(['sess ' sess])
        colormap('viridis')
        clim([0 0.27])

    case 'Spencer_2020_muscle_model'
        C = natChord_analyze('chord_distance_RDM','num_fingers',1);

        chords = C.chordID{1}(1:10);
        
        idx = [(6:10)';(1:5)'];
        chords = chords(idx);

        % calculated avg_mahalanobis:
        d_avg = 0;
        for i = 1:length(C.d_mahalanobis)
            d_avg = d_avg + d_mahalanobis/length(C.d_mahalanobis);
        end

        tmpA = avg_mahalanobis(6:10,6:10);
        tmpB = avg_mahalanobis(1:5,1:5);
        tmpC = avg_mahalanobis(1:5,6:10);
        tmpD = avg_mahalanobis(6:10,1:5);
        avg_mahalanobis(1:5,1:5) = tmpA;
        avg_mahalanobis(1:5,6:10) = tmpD;
        avg_mahalanobis(6:10,6:10) = tmpB;
        avg_mahalanobis(6:10,1:5) = tmpC;

    case 'EMG_prewhitening_matrix'
        normalize_channels = 1;
        plot_option = 1;
        subj_num = 1;
        vararginoptions(varargin,{'normalize_channels','plot_option','subj_num'})
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        data_all = dload(fullfile(project_path, 'analysis', 'natChord_all.tsv'));
        sn_unique = unique(data.sn);

        C = [];
        for sn = 1:length(sn_unique)
            % temp struct:
            C_tmp = [];
            sess = unique(data.sess(data.sn==sn_unique(sn)));
            for i = 1:length(sess)
                % calculating avg chord patterns:
                [pattern, chords] = natChord_analyze('avg_chord_patterns','subject_name',['subj' num2str(sn_unique(sn),'%.2d')],'sess',sess(i),'plot_option',0,'normalize_channels',normalize_channels);

                % make the EMG pattern residuals:
                rows = data_all.sn==sn_unique(sn) & data_all.sess==sess(i) & data_all.trialCorr==1;
                emg_patterns_all = [data_all.emg_hold_avg_e1(rows), ...
                                    data_all.emg_hold_avg_e2(rows), ...
                                    data_all.emg_hold_avg_e3(rows), ...
                                    data_all.emg_hold_avg_e4(rows), ...
                                    data_all.emg_hold_avg_e5(rows), ...
                                    data_all.emg_hold_avg_f1(rows), ...
                                    data_all.emg_hold_avg_f2(rows), ...
                                    data_all.emg_hold_avg_f3(rows), ...
                                    data_all.emg_hold_avg_f4(rows), ...
                                    data_all.emg_hold_avg_f5(rows)];
                if (normalize_channels)
                    scales = get_emg_scales(sn_unique(sn),sess(i));
                    emg_patterns_all = emg_patterns_all ./ scales;
                end

                chords_all = data_all.chordID(rows);
                pattern_residuals = zeros(size(emg_patterns_all));
                for j = 1:length(chords)
                    idx_rows = chords_all==chords(j);
                    tmp_chord_patterns = emg_patterns_all(idx_rows,:);
                    pattern_residuals(idx_rows,:) = tmp_chord_patterns - pattern(chords==chords(j),:);
                end
                
                cov_res = cov(pattern_residuals);
                C_tmp.sn = sn_unique(sn);
                C_tmp.sess = sess(i);
                C_tmp.cov_res{1} = cov_res;
                C_tmp.pattern{1} = pattern;
                C_tmp.pattern_prewhitened{1} = pattern * cov_res^-(1/2);

                C = addstruct(C,C_tmp,'row','force');
            end
            
        end 
        C.chordID = chords;
        varargout{1} = C;

        if (plot_option)
            emg_locs_names = ["e1";"e2";"e3";"e4";"e5";"f1";"f2";"f3";"f4";"f5"];
            figure('Position',[10 10 450 1500]);
            % plotting the chord EMG patterns:
            pcolor([[C.pattern{subj_num}, zeros(size(C.pattern{subj_num},1),1)] ; zeros(1,size(C.pattern{subj_num},2)+1)])
            colorbar
            % clim([0 1.5])
            
            % plot settings:
            ax = gca;
            
            set(ax,'YTick',(1:size(C.pattern{subj_num},1))+0.5)
            set(ax,'YTickLabel',chords)
            
            set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
            set(ax,'XTickLabel',emg_locs_names)
            
            set(gca,'YDir','reverse')

            figure('Position',[10 10 450 1500]);
            % plotting the chord EMG patterns:
            pcolor([[C.pattern_prewhitened{subj_num}, zeros(size(C.pattern_prewhitened{subj_num},1),1)] ; zeros(1,size(C.pattern_prewhitened{subj_num},2)+1)])
            colorbar
            % clim([0 1.5])
            
            % plot settings:
            ax = gca;
            
            set(ax,'YTick',(1:size(C.pattern_prewhitened{subj_num},1))+0.5)
            set(ax,'YTickLabel',chords)
            
            set(ax,'XTick', (1:size(emg_locs_names,1))+0.5)
            set(ax,'XTickLabel',emg_locs_names)
            
            set(gca,'YDir','reverse')
        end

    case 'emg_force_transform'
        coeff = []; % transform matrix if we want to project emg data to a certain coefficient.
        is_plot = 1;
        conference_fig = 0;
        within_sess = 0;
        vararginoptions(varargin,{'coeff','is_plot','conference_fig','within_sess'})
        
        % load data:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);
        
        [emg_rel,emg_pattern] = natChord_analyze('emg_reliability');
        [force_rel,force_pattern] = natChord_analyze('force_reliability');

        % transform the EMGs if transform matrix was given:
        if ~isempty(coeff)
            for i = 1:length(emg_pattern.pattern)
                emg_pattern.pattern{i} = emg_pattern.pattern{i} * coeff{i};
            end
        end
        
        C = [];
        for i = 1:length(emg_pattern.pattern)
            % average of in subjects:
            idx_in = setdiff(1:length(emg_pattern.pattern),i);
            if ~within_sess
                emg = 0;
                force = 0;
                for j = 1:length(idx_in)
                    emg = emg + emg_pattern.pattern{idx_in(j)}/length(idx_in);
                    force = force + force_pattern.pattern{idx_in(j)}/length(idx_in);
                end
            else
                emg = emg_pattern.pattern{i};
                force = force_pattern.pattern{i};
            end
            directions = make_design_matrix(chords,'additive');
            
            % EMG to force linear mixture matrix:
            W = (emg' * emg)^-1 * emg' * force;

            % force to EMG linear mixture matrix:
            A = (force' * force)^-1 * force' * emg;

            % EMG to behavior linear mixture matrix:
            W_b = (emg' * emg)^-1 * emg' * directions;
            
            % behavior to EMG linear mixture matrix:
            A_b = (directions' * directions)^-1 * directions' * emg;

            % calculate the R^2 of reconstructed values:
            force_pred = emg*W;
            R2_m2f = mean(1 - sum((force_pred-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            R_m2f = corr2(force_pred,force_pattern.pattern{i});
            
            directions_pred = emg*W_b;
            R2_m2d = mean(1 - sum((directions_pred-directions).^2,1) ./ sum(directions.^2,1));
            R_m2d = corr2(directions_pred,directions);
            
            emg_pred = force*A;
            R2_f2m = mean(1 - sum((emg_pred-emg_pattern.pattern{i}).^2,1) ./ sum(emg_pattern.pattern{i}.^2,1));
            R_f2m = corr2(emg_pred,emg_pattern.pattern{i});
            
            emg_pred = directions*A_b;
            R2_d2m = mean(1 - sum((emg_pred-emg_pattern.pattern{i}).^2,1) ./ sum(emg_pattern.pattern{i}.^2,1));
            R_d2m = corr2(emg_pred,emg_pattern.pattern{i});
            
            tmp.R2_m2f = R2_m2f;
            tmp.R2_m2d = R2_m2d;
            tmp.R2_f2m = R2_f2m;
            tmp.R2_d2m = R2_d2m;
            tmp.R_m2f = R_m2f;
            tmp.R_m2d = R_m2d;
            tmp.R_f2m = R_f2m;
            tmp.R_d2m = R_d2m;

            % cross validation with 1-chord out:
            F2M_crossval = zeros(size(emg));
            D2M_crossval = zeros(size(emg));
            F_crossval = zeros(size(force));
            D_crossval = zeros(size(directions));
            for ch = 1:length(chords)
                % remove one chord from the train matrices:
                emg_reduced = emg(setdiff(1:length(chords),ch),:);
                force_reduced = force(setdiff(1:length(chords),ch),:);
                directions_reduced = directions(setdiff(1:length(chords),ch),:);
            
                % EMG to force linear mixture matrix:
                W = (emg_reduced' * emg_reduced)^-1 * emg_reduced' * force_reduced;
                % force to EMG linear mixture matrix:
                A = (force_reduced' * force_reduced)^-1 * force_reduced' * emg_reduced;
                % EMG to directions linear mixture matrix:
                W_b = (emg_reduced' * emg_reduced)^-1 * emg_reduced' * directions_reduced;
                % directions to EMG linear mixture matrix:
                A_b = (directions_reduced' * directions_reduced)^-1 * directions_reduced' * emg_reduced;

                % add predictions to the corresponding matrices:
                force_pred = emg(ch,:)*W;
                F_crossval(ch,:) = force_pred;

                directions_pred = emg(ch,:)*W_b;
                D_crossval(ch,:) = directions_pred;

                emg_pred = force(ch,:)*A;
                F2M_crossval(ch,:) = emg_pred;

                emg_pred = directions(ch,:)*A_b;
                D2M_crossval(ch,:) = emg_pred;
                
            end
            tmp.R2_m2f_crossval = mean(1 - sum((F_crossval-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            tmp.R2_m2d_crossval = mean(1 - sum((D_crossval-directions).^2,1) ./ sum(directions.^2,1));
            tmp.R2_f2m_crossval = mean(1 - sum((F2M_crossval-emg_pattern.pattern{i}).^2,1) ./ sum(emg_pattern.pattern{i}.^2,1));
            tmp.R2_d2m_crossval = mean(1 - sum((D2M_crossval-emg_pattern.pattern{i}).^2,1) ./ sum(emg_pattern.pattern{i}.^2,1));
            tmp.R_m2f_crossval = corr2(F_crossval,force_pattern.pattern{i});
            tmp.R_m2d_crossval = corr2(D_crossval,directions);
            tmp.R_f2m_crossval = corr2(F2M_crossval,emg_pattern.pattern{i});
            tmp.R_d2m_crossval = corr2(D2M_crossval,emg_pattern.pattern{i});
            
            C = addstruct(C,tmp,'row','force');
        end
        varargout{1} = C;
        
        % significant tests:
        % [t,p] = ttest(C.R2_force2emg_crossval,C.R2_emg2force_crossval,1,'paired');
        % fprintf('test if force2emg > emg2force: %.4f\n',p)
        % [t,p] = ttest(C.R2_description2emg_crossval,C.R2_emg2description_crossval,1,'paired');
        % fprintf('test if direction2emg > emg2direction: %.4f\n',p)
        
        if is_plot
            if ~ conference_fig
                % PLOT 1:
                figure;
                ax1 = axes('Units', 'centimeters', 'Position', [2 2 4.8 5],'Box','off');
                ax1.PositionConstraint = "innerposition";
                axes(ax1);
                x = [[ones(length(C.R2_f2m),1) ; 2*ones(length(C.R2_m2f),1)] , [4*ones(length(C.R2_f2m_crossval),1) ; 5*ones(length(C.R2_m2f_crossval),1)]];
                y = [[C.R2_f2m ; C.R2_f2m_crossval] , [C.R2_m2f ; C.R2_m2f_crossval]];
                [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',1,'gapwidth',[1,0,0,0]); 
                lim_width = 0.4;
                drawline(mean(emg_rel.r2),'dir','horz','lim',[x_coord(1)-lim_width,x_coord(1)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(force_rel.r2),'dir','horz','lim',[x_coord(2)-lim_width,x_coord(2)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(emg_rel.r2),'dir','horz','lim',[x_coord(3)-lim_width,x_coord(3)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(force_rel.r2),'dir','horz','lim',[x_coord(4)-lim_width,x_coord(4)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
        
                ylim([0.5,1])
                xlim([0,6])
                h = gca;
                h.XTickLabels = {'Force to EMG','EMG to Force','Force to EMG','EMG to Force'};
                h.XAxis.FontSize = my_font.xlabel;
                h.YAxis.FontSize = my_font.tick_label;
                h.LineWidth = 2;
                ylabel('$\mathbf{R^2}$','interpreter','LaTex','FontSize',my_font.ylabel)
                fontname("Arial")
        
                % PLOT 2:
                figure;
                ax1 = axes('Units', 'centimeters', 'Position', [2 2 4.8 5],'Box','off');
                ax1.PositionConstraint = "innerposition";
                axes(ax1);
                x = [[ones(length(C.R_f2m),1) ; 2*ones(length(C.R_m2f),1)] , [4*ones(length(C.R_f2m_crossval),1) ; 5*ones(length(C.R_m2f_crossval),1)]];
                y = [[C.R_f2m ; C.R_f2m_crossval] , [C.R_m2f ; C.R_m2f_crossval]];
                [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',1,'gapwidth',[1,0,0,0]); 
                lim_width = 0.4;
                drawline(mean(emg_rel.r),'dir','horz','lim',[x_coord(1)-lim_width,x_coord(1)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(force_rel.r),'dir','horz','lim',[x_coord(2)-lim_width,x_coord(2)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(emg_rel.r),'dir','horz','lim',[x_coord(3)-lim_width,x_coord(3)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                drawline(mean(force_rel.r),'dir','horz','lim',[x_coord(4)-lim_width,x_coord(4)+lim_width],'linewidth',2,'color',[0.85 0.85 0.85],'linestyle',':')
                
                ylim([0.5,1])
                xlim([0,6])
                h = gca;
                h.XTickLabels = {'Force to EMG','EMG to Force','Force to EMG','EMG to Force'};
                h.XAxis.FontSize = my_font.xlabel;
                h.YAxis.FontSize = my_font.tick_label;
                h.LineWidth = 2;
                ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.ylabel)
                fontname("Arial")
            else
                fig = figure('Units','centimeters', 'Position',[15 15 12 15]);
                x = [ones(length(C.R_f2m),1) ; 2*ones(length(C.R_m2f),1)];
                y = [C.R_f2m ; C.R_m2f];
                [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_line_width,'gapwidth',[0.5,0,0,0],'barwidth',1); 
                lim_width = 0.5;
                drawline(mean(emg_rel.r),'dir','horz','lim',[x_coord(1)-lim_width,x_coord(1)+lim_width],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
                drawline(mean(force_rel.r),'dir','horz','lim',[x_coord(2)-lim_width,x_coord(2)+lim_width],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
                % drawline(mean(emg_rel.r),'dir','horz','lim',[x_coord(3)-lim_width,x_coord(3)+lim_width],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
                % drawline(mean(force_rel.r),'dir','horz','lim',[x_coord(4)-lim_width,x_coord(4)+lim_width],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
                % 
                ylim([0,1])
                % xlim([0,6])
                h = gca;
                h.XTickLabels = {'Force to EMG','EMG to Force','Force to EMG','EMG to Force'};
                h.XAxis.FontSize = my_font.conf_label;
                h.YAxis.FontSize = my_font.conf_tick_label;
                h.LineWidth = conf.axis_width;
                ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
                fontname("Arial")
            end
        end

    case 'forward_selection'
        % handling input arguments:
        alpha = 0.05;
        measure = 'MD';
        sess = [3,4];
        models = {'n_fing','n_fing+additive','n_fing+2fing_adj','n_fing+2fing','n_fing+nSphere_avg','n_fing+magnitude_avg','n_fing+emg_additive_avg','n_fing+emg_2channel_avg','n_fing+force_avg','n_fing+force_2fing'};
        % models = {'n_fing+emg_additive_avg','n_fing+emg_2channel_avg','n_fing+magnitude_avg','n_fing+nSphere_avg'};
        % models = {'n_fing+additive','n_fing+2fing_adj','n_fing+2fing'};
        vararginoptions(varargin,{'alpha','measure','models','sess'})
        base_models = models;
        
        % loading data:
        chords_natChord = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords_natChord = chords_natChord.chordID(chords_natChord.sn==1 & chords_natChord.sess==1);
        data = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        data = getrow(data,ismember(data.chordID,chords_natChord));
        chords = data.chordID(data.sn==1 & data.sess==1);
        subj = unique(data.sn);
        
        % getting the values of measure:
        values_tmp = eval(['data.' measure]);

        % getting the average of sessions for every subj:
        values = zeros(length(chords),length(subj));
        for i = 1:length(subj)
            % avg with considering nan values since subjects might have
            % missed all 5 repetitions in one session:
            values(:,i) = mean([values_tmp(data.sess==sess(1) & data.sn==subj(i)),values_tmp(data.sess==sess(2) & data.sn==subj(i))],2,'omitmissing');
        end
        
        % selection steps:
        steps = length(models);
        winning_model = '';
        best_r = zeros(length(subj),1);
        C = [];
        for i = 1:steps
            for j = 1:length(models)
                r = zeros(length(subj),1);
                % loop on subjects and regression with leave-one-out:
                for sn = 1:length(subj)
                    % values of 'in' subjects, Nx1 vector:
                    y_train = values(:,setdiff(1:length(subj),sn));
                    y_train = mean(y_train,2);
        
                    % avg of 'out' subject:
                    y_test = values(:,sn);

                    % getting design matrix for model:
                    X = make_design_matrix(chords,models{j});
    
                    % training the model:
                    % [B,STATS] = linregress(y_train,X,'intercept',0);
                    [B,~] = svd_linregress(y_train,X);
    
                    % testing the model:
                    X_test = make_design_matrix(chords,models{j});
                    y_pred = X_test * B;
                    r(sn) = corr(y_pred,y_test);
                end
                tmp.step = i;
                tmp.model = models(j);
                tmp.r = {r};
                tmp.r_avg = mean(r);
                [~,tmp.pval] = ttest(r,best_r,1,'paired');
                tmp.significant = tmp.pval < alpha;

                % initialize this variable for later steps:
                tmp.win = 0;

                % save values in struct:
                C = addstruct(C,tmp,'row','force');
            end

            % competition between models:
            r_avg = C.r_avg(C.step==i);
            p_val = C.pval(C.step==i);

            % find the significant improvements:
            p_val = p_val < alpha;
            
            % if there was at least one significantly better model:
            if sum(p_val)~=0
                % remove the non-significant models from the comptetition:
                r_avg(p_val==0) = 0;
            end
            
            % find the best (significant) model:
            [~,idx] = sort(r_avg);
            
            % conclude the winner and give it a gold medal:
            winning_model = models{idx(end)};
            C.win(C.step==i & strcmp(C.model,winning_model)) = 1;
            % set the competition values for the next step:
            best_r = C.r{C.step==i & strcmp(C.model,winning_model)};

            % make models for the next step:
            models = base_models;
            split_names = strsplit(winning_model,'+');
            for j = 1:length(split_names)
                models(strcmp(models,split_names{j})) = [];
            end
            prefix = [winning_model , '+'];
            models = cellfun(@(x) [prefix x], models, 'UniformOutput', false);
        end
        
        % fing the significant winner:
        significant_steps = zeros(length(unique(C.step)),1);
        for i = 1:length(unique(C.step))
            significant = C.significant(C.step==i);
            if sum(significant)
                significant_steps(i) = 1;
            end
        end
        idx = find(significant_steps);
        idx = idx(end);
        significant_winner = C.model{C.step==idx & C.win==1};

        fprintf('\nThe significant winner of the forward selection is:\n%s\n',significant_winner)
        fprintf('\nThe r_avg winner of forward selection is:\n%s\n',winning_model);
        varargout{1} = C;

    case 'backward_selection'
        % handling input arguments:
        alpha = 0.05;
        measure = 'MD';
        sess = [3,4];
        models = {'additive','2fing_adj','2fing','nSphere_avg','magnitude_avg','emg_additive_avg','emg_2channel_avg','force_avg','force_2fing'};
        vararginoptions(varargin,{'alpha','measure','models','sess'})
        base_models = models;
        
        % loading data:
        chords_natChord = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords_natChord = chords_natChord.chordID(chords_natChord.sn==1 & chords_natChord.sess==1);
        data = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        data = getrow(data,ismember(data.chordID,chords_natChord));
        chords = data.chordID(data.sn==1 & data.sess==1);
        subj = unique(data.sn);
        
        % getting the values of measure:
        values_tmp = eval(['data.' measure]);
        
        % getting the average of sessions for every subj:
        values = zeros(length(chords),length(subj));
        for i = 1:length(subj)
            % avg with considering nan values since subjects might have
            % missed all 5 repetitions in one session:
            values(:,i) = mean([values_tmp(data.sess==sess(1) & data.sn==subj(i)),values_tmp(data.sess==sess(2) & data.sn==subj(i))],2,'omitmissing');
        end
        
        % full model:
        for i = 1:length(base_models)
            full_model = strjoin(base_models, '+');
        end
        
        best_r = zeros(length(subj),1);
        % loop on subjects and regression with leave-one-out:
        for sn = 1:length(subj)
            % values of 'in' subjects, Nx1 vector:
            y_train = values(:,setdiff(1:length(subj),sn));
            y_train = mean(y_train,2);

            % avg of 'out' subject:
            y_test = values(:,sn);

            % getting design matrix for model:
            X = make_design_matrix(chords,full_model);

            % training the model:
            % [B,STATS] = linregress(y_train,X,'intercept',0);
            [B,~] = svd_linregress(y_train,X);

            % testing the model:
            X_test = make_design_matrix(chords,full_model);
            y_pred = X_test * B;
            best_r(sn) = corr(y_pred,y_test);
        end
        
        % selection steps:
        steps = length(models)-1;

        % define models for the first step:
        winning_model = full_model;
        
        C = [];
        for i = 1:steps
            % define models to be removed for the current step:
            reduce_models = strsplit(winning_model,'+');
            
            for j = 1:length(reduce_models)
                model = reduce_models;
                model(strcmp(model,reduce_models{j})) = [];
                model = strjoin(model,'+');

                r = zeros(length(subj),1);
                % loop on subjects and regression with leave-one-out:
                for sn = 1:length(subj)
                    % values of 'in' subjects, Nx1 vector:
                    y_train = values(:,setdiff(1:length(subj),sn));
                    y_train = mean(y_train,2);
        
                    % avg of 'out' subject:
                    y_test = values(:,sn);

                    % getting design matrix for model:
                    X = make_design_matrix(chords,model);
    
                    % training the model:
                    % [B,STATS] = linregress(y_train,X,'intercept',0);
                    [B,~] = svd_linregress(y_train,X);
    
                    % testing the model:
                    X_test = make_design_matrix(chords,model);
                    y_pred = X_test * B;
                    r(sn) = corr(y_pred,y_test);
                end
                tmp.step = i;
                tmp.model = {model};
                tmp.r = {r};
                tmp.r_avg = mean(r);
                [~,tmp.pval] = ttest(best_r,r,1,'paired');
                tmp.significantly_worse = tmp.pval < alpha; % this is 1 if the reduced model is significantly worse than the starting model
                tmp.win = 0;

                if tmp.significantly_worse == 0
                    tmp.fail = 1;   % the reduced model fails if it does not significantly reduce performance
                end
                
                % save values in struct:
                C = addstruct(C,tmp,'row','force');
            end

            % conclude the step and choose a winning model for the next step:
            r_avg = C.r_avg(C.step==i);
            fail = C.fail(C.step==i);
            
            % if there was at least one significantly worse model:
            if sum(fail) ~= length(fail)
                % remove the good models from the competition
                r_avg(fail==0) = 0;
            end
            
            % sort the reduction of performance after removing each model:
            [~,idx] = sort(mean(best_r)-r_avg);

            % conclude the not-loser(winner?) and give it a gold medal:
            tmp_models = C.model(C.step==i);
            winning_model = tmp_models{idx(1)};
            C.win(C.step==i & strcmp(C.model,winning_model)) = 1;
            
            % set the competition values for the next step:
            best_r = C.r{C.step==i & strcmp(C.model,winning_model)};
        end
        
        % find the significant winner:
        % significant_winner = '';
        % significant_steps = zeros(length(unique(C.step)),1);
        % for i = 1:length(unique(C.step))
        %     significant = C.significant(C.step==i);
        %     if sum(significant)
        %         significant_steps(i) = 1;
        %     end
        % end
        % idx = find(significant_steps);
        % idx = idx(end);
        % significant_winner = C.model{C.step==idx & C.win==1};
        % 
        % fprintf('\nThe significant winner of the forward selection is:\n%s\n',significant_winner)
        % fprintf('\nThe r_avg winner of forward selection is:\n%s\n',winning_model);
        varargout{1} = C;
        
    case 'emg_reliability'
        is_plot = 0;
        vararginoptions(varargin,'is_plot');
        data = dload(fullfile(project_path,'analysis','natChord_analysis.tsv'));
        subj = unique(data.sn);
        
        % save all emg patterns:
        emg_pattern = [];
        for i = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(i)));
            for j = 1:length(sess)
                tmp = [];
                row = data.sn==subj(i) & data.sess==sess(j);
                tmp_pattern = [data.emg_e1(row),data.emg_e2(row),data.emg_e3(row),data.emg_e4(row),data.emg_e5(row),...
                               data.emg_f1(row),data.emg_f2(row),data.emg_f3(row),data.emg_f4(row),data.emg_f5(row)];
                tmp.sn = subj(i);
                tmp.sess = sess(j);
                tmp.pattern = {tmp_pattern};
                emg_pattern = addstruct(emg_pattern,tmp,'row','force');
            end
        end

        % estimate leave-one-out reliability
        C = [];
        for i = 1:length(emg_pattern.pattern)
            idx_in = setdiff(1:length(emg_pattern.pattern),i);
            avg = 0;
            for j = 1:length(idx_in)
                avg = avg + emg_pattern.pattern{idx_in(j)}/length(idx_in);
            end
            
            tmp = [];
            % corr of leave subject with average of included subjects:
            tmp.r = corr2(avg,emg_pattern.pattern{i});
            % r2 of out subejct with  average of included subjects:
            SSE = sum((emg_pattern.pattern{i} - avg).^2,1);
            SST = sum(avg.^2,1);
            tmp.r2 = mean(1 - SSE./SST);
            C = addstruct(C,tmp,'row','force');
        end
        varargout{1} = C;
        varargout{2} = emg_pattern;

        % PLOTS:
        if is_plot
            avg_pattern = mean(cat(3, emg_pattern.pattern{:}), 3)';
            % swap flexors and extensor electrodes location in matrix:
            swap = avg_pattern(6:10,:);
            avg_pattern(6:10,:) = avg_pattern(1:5,:);
            avg_pattern(1:5,:) = swap;

            % swap single finger flexion and extensions in matrix:
            swap = avg_pattern(:,6:10);
            avg_pattern(:,6:10) = avg_pattern(:,1:5);
            avg_pattern(:,1:5) = swap;

            fig = figure('Units','centimeters', 'Position',[15 15 20 14]);
            pattern_1f = avg_pattern(:,[16,35,62,68]);
            p = pcolor([[pattern_1f, zeros(size(pattern_1f,1),1)] ; zeros(1,size(pattern_1f,2)+1)]);
            colormap('viridis')
            clim([0 0.27])
            colorbar
            % plot settings:
            set(p, 'EdgeColor', [0.2,0.2,0.2]);
            set(gca,'YDir','reverse')
            set(gca,'XTick', (1:size(pattern_1f,1))+0.5)
            set(gca,'YTick', (1:size(pattern_1f,2))+0.5)
            % axis square
            % 
            % xlabel('days','FontSize',my_font.conf_label)
            % ylabel(['gooz'],'FontSize',my_font.conf_label)
            % % ylabel([measure],'FontSize',my_font.tick_label)
            % h = gca;
            % h.YTick = linspace(h.YTick(1),h.YTick(end),3);
            % h.XAxis.FontSize = my_font.conf_tick_label;
            % h.YAxis.FontSize = my_font.conf_tick_label;
            % h.LineWidth = conf.axis_width;
            fontname("Arial")
        end

    case 'force_reliability'
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        subj = unique(data.sn);
        
        % save all emg patterns:
        force_pattern = [];
        for i = 1:length(subj)
            sess = unique(data.sess(data.sn==subj(i)));
            for j = 1:length(sess)
                tmp = [];
                row = data.sn==subj(i) & data.sess==sess(j);
                diff_force = [data.diff_force_f1(row),...
                              data.diff_force_f2(row),...
                              data.diff_force_f3(row),...
                              data.diff_force_f4(row),...
                              data.diff_force_f5(row)];
                force = [diff_force,diff_force];
                for col = 1:5
                    force(force(:,col)<0,col) = 0;
                    force(force(:,col+5)>0,col+5) = 0;
                end
                force(:,6:10) = force(:,6:10)*-1;


                tmp.sn = subj(i);
                tmp.sess = sess(j);
                tmp.pattern = {force};
                force_pattern = addstruct(force_pattern,tmp,'row','force');
            end
        end

        % estimate leave-one-out reliability
        C = [];
        for i = 1:length(force_pattern.pattern)
            idx_in = setdiff(1:length(force_pattern.pattern),i);
            avg = 0;
            for j = 1:length(idx_in)
                avg = avg + force_pattern.pattern{idx_in(j)}/length(idx_in);
            end
            
            tmp = [];
            % corr of leave subject with average of included subjects:
            tmp.r = corr2(avg,force_pattern.pattern{i});
            % r2 of out subejct with  average of included subjects:
            SSE = sum((force_pattern.pattern{i} - avg).^2,1);
            SST = sum(avg.^2,1);
            tmp.r2 = mean(1 - SSE./SST);
            C = addstruct(C,tmp,'row','force');
        end
        varargout{1} = C;
        varargout{2} = force_pattern;
        
    case 'emg_behaviour_model_comparison'
        conference_fig = 0;
        vararginoptions(varargin,{'conference_fig'})
        model_names = {'n_fing+additive','n_fing+emg_additive_avg','n_fing+additive+2fing','n_fing+emg_additive_avg+emg_2channel_avg'};
        [C,stats] = natChord_analyze('model_testing_all_efc1','measure','MD','model_names',model_names,'is_plot',0);
        varargout{1} = C;
        varargout{2} = stats;

        C_MD = natChord_analyze('model_testing_all_efc1','measure','MD','model_names',{'n_fing'},'is_plot',0);
        r_finger_count = mean(C_MD.r);
        
        % PLOT:
        if ~conference_fig
            figure;
            ax1 = axes('Units', 'centimeters', 'Position', [2 2 4.8 5],'Box','off');
            ax1.PositionConstraint = "innerposition";
            axes(ax1);
            drawline(C.noise_ceil(1),'dir','horz','lim',[0,6],'linewidth',2,'color',[0.7 0.7 0.7],'linestyle',':'); hold on;

            N = length(C.r(strcmp(C.model,model_names{1})));
            x = [[ones(N,1) ; 2*ones(N,1)] , [4*ones(N,1) ; 5*ones(N,1)]];
            y = [[C.r(strcmp(C.model,model_names{1})) ; C.r(strcmp(C.model,model_names{3}))] , [C.r(strcmp(C.model,model_names{2})) ; C.r(strcmp(C.model,model_names{4}))]];
            [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',1,'gapwidth',[1,0,0,0]); 
            
            ylim([r_finger_count,1])
            xlim([0,6])
            h = gca;
            h.XTick = [1.5 , 4.5];
            h.YTick = 0.5:0.1:1;
            h.XTickLabels = {'Additive','Additive+Interaction'};
            h.XAxis.FontSize = my_font.xlabel;
            h.YAxis.FontSize = my_font.tick_label;
            ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.ylabel)
            fontname("Arial")
        else
            fig = figure('Units','centimeters', 'Position',[15 15 23 24]);
            drawline(C.noise_ceil(1),'dir','horz','lim',[0,6],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':'); hold on;
            N = length(C.r(strcmp(C.model,model_names{1})));
            x = [[ones(N,1) ; 2*ones(N,1)] , [4*ones(N,1) ; 5*ones(N,1)]];
            y = [[C.r(strcmp(C.model,model_names{1})) ; C.r(strcmp(C.model,model_names{3}))] , [C.r(strcmp(C.model,model_names{2})) ; C.r(strcmp(C.model,model_names{4}))]];
            [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',3,'gapwidth',[1,0,0,0],'barwidth',1); 
            ylim([r_finger_count,0.89])
            xlim([0,6])
            h = gca;
            h.XTick = [1.5 , 4.5];
            h.YTick = [0.735,round(r_finger_count/2,2)+round(C.noise_ceil(1)/2,2),round(C.noise_ceil(1),2)];
            h.XTickLabels = {'Additive','Additive+Interaction'};
            h.XAxis.FontSize = my_font.conf_label;
            h.YAxis.FontSize = my_font.conf_tick_label;
            h.LineWidth = conf.axis_width;
            ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
            fontname("Arial")
        end

    case 'PCA_emg'
        sampling_option = 'whole_sampled';
        num_dim = 4;
        vararginoptions(varargin,{'sampling_option','num_dim'});
        
        % loading data:
        [emg_rel,emg_pattern] = natChord_analyze('emg_reliability');
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);
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
            
            sess = unique(data.sess(data.sn == subjects(sn)));
            for j = 1:length(sess)
                % scaling factors:
                scales = get_emg_scales(subjects(sn),sess(j));
                
                emg_dist_sess = getrow(emg_dist,emg_dist.sess==sess(j));
                % normalizing the natural EMGs:
                for i = 1:length(emg_dist_sess.dist)
                    emg_dist_sess.dist{i} = (emg_dist_sess.dist{i} ./ scales);
                end
                
                % getting avg mean deviation of chords across sessions:
                MD = zeros(length(chords),1);
                for k = 1:length(chords)
                    row = data.sn==subjects(sn) & data.chordID==chords(k);
                    MD(k) = data.MD_efc(row & data.sess==sess(j));
                end
                
                % PCA of natural and chord EMG:
                natural_emg = vertcat(emg_dist_sess.dist{:});
                [coeff_nat, score_nat, eigenvalues_nat, ~, explained_nat] = pca(natural_emg);
                [coeff_chord, score_chord, eigenvalues_chord, ~, explained_chord] = pca(emg_pattern.pattern{emg_pattern.sn==subjects(sn) & emg_pattern.sess==sess(j)});
                tmp.sn = subjects(sn);
                tmp.sess = sess(j);
                tmp.explained_nat = explained_nat';
                tmp.explained_chord = explained_chord';
                tmp.coeff_nat = {coeff_nat};
                tmp.coeff_chord = {coeff_chord};
                tmp.chord_pattern = {emg_pattern.pattern{emg_pattern.sn==subjects(sn) & emg_pattern.sess==sess(j)}};
                tmp.natural_emg = {natural_emg};
                C = addstruct(C,tmp,'row','force');
            end
        end
        varargout{1} = C;

        % % force explained by projected EMG to natural and Chord:
        % C0 = natChord_analyze('emg_force_transform','is_plot',0);
        % C_nat = natChord_analyze('emg_force_transform','is_plot',0,'coeff',C.coeff_nat);
        % C_chord = natChord_analyze('emg_force_transform','is_plot',0,'coeff',C.coeff_chord);
        % % PLOT:
        % fig = figure('Units','centimeters', 'Position',[15 15 25 25]);
        % x = [ones(length(C_nat.R_m2f),1) ; 2*ones(length(C_chord.R_m2f),1)];
        % y = [C_nat.R_m2f ; C_chord.R_m2f];
        % [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_width); 
        % drawline(mean(C0.R_m2f),'dir','horz','lim',[0,3.5],'linewidth',4,'color',[0.85 0.85 0.85],'linestyle',':')
        % 
        % ylim([0,1])
        % xlim([0,3.5])
        % h = gca;
        % h.XTickLabels = {'Project to Natural','Project to Chord'};
        % h.XAxis.FontSize = my_font.conf_tick_label;
        % h.YAxis.FontSize = my_font.conf_tick_label;
        % h.LineWidth = 2;
        % ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
        % fontname("Arial")
        % 
        % % projecting the EMG to different number of coeffs and explaining force:
        % for i = 1:length(C.coeff_nat)
        %     coeff1{i} = [C.coeff_nat{i}(:,1:num_dim)];
        %     coeff2{i} = [C.coeff_nat{i}(:,num_dim+1:end)];
        % end
        % C1 = natChord_analyze('emg_force_transform','coeff',coeff1,'is_plot',0);
        % C2 = natChord_analyze('emg_force_transform','coeff',coeff2,'is_plot',0);
        % % PLOT:
        % fig = figure('Units','centimeters', 'Position',[15 15 25 25]);
        % x = [ones(length(C1.R_m2f),1) ; 2*ones(length(C2.R_m2f),1)];
        % y = [C1.R_m2f ; C2.R_m2f];
        % [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_width); 
        % drawline(mean(C_nat.R_m2f),'dir','horz','lim',[0,3.5],'linewidth',4,'color',[0.85 0.85 0.85],'linestyle',':')
        % 
        % y = C.explained_nat';
        % eig_sum = sum(y,1);
        % y = y./eig_sum;
        % y = cumsum(y,1)*100;
        % y = mean(y,2);
        % 
        % ylim([0,1])
        % xlim([0,3.5])
        % h = gca;
        % h.XTickLabels = {num2str(y(num_dim),'%.2f'),num2str(100-y(num_dim),'%.2f')};
        % h.XAxis.FontSize = my_font.conf_tick_label;
        % h.YAxis.FontSize = my_font.conf_tick_label;
        % h.LineWidth = 2;
        % ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
        % fontname("Arial")
        % 
        % % PLOT - var explained by dimension:
        % fig = figure('Units','centimeters', 'Position',[15 15 15 15]);
        % drawline(100,'dir','horz','lim',[0,11],'linewidth',4,'color',[0.8 0.8 0.8],'linestyle','-'); hold on;
        % 
        % x = repmat((1:10)',size(C.explained_nat,1),1);
        % y = C.explained_nat';
        % eig_sum = sum(y,1);
        % y = y./eig_sum;
        % y = cumsum(y,1)*100;
        % y = y(:);
        % lineplot(x,y,'markertype','o','markersize',5,'markerfill',colors_red(2,:),'markercolor',colors_red(2,:),'linecolor',colors_red(2,:),'linewidth',2,'errorcolor',colors_red(2,:)); hold on;
        % 
        % x = repmat((1:10)',size(C.explained_chord,1),1);        % y = C.explained_chord';
        % eig_sum = sum(y,1);
        % y = y./eig_sum;
        % y = cumsum(y,1)*100;
        % y = y(:);
        % lineplot(x,y,'markersize',5,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',2,'errorcolor',colors_blue(5,:));
        % lgd = legend({'','natural','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','chord'});
        % xlim([0,11])
        % ylim([0,105])
        % ylabel('Variance Explained (%)')
        % xlabel('Number of Dimensions')

    case 'PCA_emg_crippled'
        num_dim = 4;
        project_to_chord = 0;
        vararginoptions(varargin,{'num_dim','project_to_chord'})
        % load data:
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);
        
        % get PCA coeffs:
        C_PCA = natChord_analyze('PCA_emg');
        
        % variance explained by PCs:
        if ~project_to_chord
            var1 = mean(sum(C_PCA.explained_nat(:,1:num_dim),2));
            var2 = mean(sum(C_PCA.explained_nat(:,num_dim+1:end),2));
        else
            var1 = mean(sum(C_PCA.explained_chord(:,1:num_dim),2));
            var2 = mean(sum(C_PCA.explained_chord(:,num_dim+1:end),2));
        end
        
        [emg_rel,emg_pattern] = natChord_analyze('emg_reliability');
        [force_rel,force_pattern] = natChord_analyze('force_reliability');
        
        C = [];
        for i = 1:length(emg_pattern.pattern)
            emg_nat = C_PCA.natural_emg{i};
            emg_nat = emg_nat - mean(emg_nat,1);
            emg = emg_pattern.pattern{i};
            emg = emg - mean(emg,1);
            if ~project_to_chord
                emg_proj = emg * C_PCA.coeff_nat{i}; % project EMG to natural PCs
                emg_proj = emg_proj-mean(emg_proj,1); 
                emg_crippled01 = emg * C_PCA.coeff_nat{i}(:,1:num_dim);
                emg_crippled01 = emg_crippled01 - mean(emg_crippled01,1);
                emg_crippled02 = emg * C_PCA.coeff_nat{i}(:,num_dim+1:end);
                emg_crippled02 = emg_crippled02 - mean(emg_crippled02,1);
            else
                emg_proj = emg * C_PCA.coeff_chord{i}; % project EMG to natural PCs
                emg_proj = emg_proj-mean(emg_proj,1); 
                emg_crippled01 = emg * C_PCA.coeff_chord{i}(:,1:num_dim);
                emg_crippled01 = emg_crippled01 - mean(emg_crippled01,1);
                emg_crippled02 = emg * C_PCA.coeff_chord{i}(:,num_dim+1:end);
                emg_crippled02 = emg_crippled02 - mean(emg_crippled02,1);s
            end
            
            force = force_pattern.pattern{i};
            force = force - mean(force,1);

            % measure the variance of EMG after projection on each eigen vector:
            for j = 1:10
                M = emg'*emg;
                M_nat = emg_nat' * emg_nat;

                total_var_chord = trace(M);
                total_var_nat = trace(M_nat);

                total_var_chord2nat = trace(C_PCA.coeff_nat{i}' * M * C_PCA.coeff_nat{i});
                total_var_chord2chord = trace(C_PCA.coeff_chord{i}' * M * C_PCA.coeff_chord{i});

                total_var_nat2nat = trace(C_PCA.coeff_nat{i}' * M_nat * C_PCA.coeff_nat{i});
                total_var_nat2chord = trace(C_PCA.coeff_chord{i}' * M_nat * C_PCA.coeff_chord{i});
                
                tmp.var_chord2nat(1,j) = trace(C_PCA.coeff_nat{i}(:,j)' * M * C_PCA.coeff_nat{i}(:,j));
                tmp.var_chord2chord(1,j) = trace(C_PCA.coeff_chord{i}(:,j)' * M * C_PCA.coeff_chord{i}(:,j));
                
                tmp.var_nat2nat(1,j) = trace(C_PCA.coeff_nat{i}(:,j)' * M_nat * C_PCA.coeff_nat{i}(:,j));
                tmp.var_nat2chord(1,j) = trace(C_PCA.coeff_chord{i}(:,j)' * M_nat * C_PCA.coeff_chord{i}(:,j));
            end
            tmp.total_var_chord = total_var_chord;
            tmp.total_var_nat = total_var_nat;
            tmp.total_var_chord2nat = total_var_chord2nat;
            tmp.total_var_chord2chord = total_var_chord2chord;
            tmp.total_var_nat2nat = total_var_nat2nat;
            tmp.total_var_nat2chord = total_var_nat2chord;
            
            % EMG to force linear mixture matrix:
            W = (emg' * emg)^-1 * emg' * force;
            W_proj = (emg_proj' * emg_proj)^-1 * emg_proj' * force;
            W_crippled01 = (emg_crippled01' * emg_crippled01)^-1 * emg_crippled01' * force;
            W_crippled02 = (emg_crippled02' * emg_crippled02)^-1 * emg_crippled02' * force;
            % calculate the R^2 of reconstructed force:
            force_pred = emg*W;
            force_pred_proj = emg_proj * W_proj;
            force_pred_crippled01 = emg_crippled01 * W_crippled01;
            force_pred_crippled02 = emg_crippled02 * W_crippled02;

            R2_m2f = mean(1 - sum((force_pred-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            R_m2f = corr2(force_pred,force_pattern.pattern{i});

            R2_m2f_proj = mean(1 - sum((force_pred_proj-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            R_m2f_proj = corr2(force_pred_proj,force_pattern.pattern{i});

            R2_m2f_crippled01 = mean(1 - sum((force_pred_crippled01-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            R_m2f_crippled01 = corr2(force_pred_crippled01,force_pattern.pattern{i});

            R2_m2f_crippled02 = mean(1 - sum((force_pred_crippled02-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            R_m2f_crippled02 = corr2(force_pred_crippled02,force_pattern.pattern{i});
            
            tmp.R2_m2f = R2_m2f;
            tmp.R_m2f = R_m2f;
            tmp.R2_m2f_proj = R2_m2f_proj;
            tmp.R_m2f_proj = R_m2f_proj;
            tmp.R2_m2f_crippled01 = R2_m2f_crippled01;
            tmp.R_m2f_crippled01 = R_m2f_crippled01;
            tmp.R2_m2f_crippled02 = R2_m2f_crippled02;
            tmp.R_m2f_crippled02 = R_m2f_crippled02;

            % cross validation with 1-chord out:
            F_crossval = zeros(size(force));
            F_crossval_proj = zeros(size(force));
            F_crossval_crippled01 = zeros(size(force));
            F_crossval_crippled02 = zeros(size(force));
            for ch = 1:length(chords)
                % remove one chord from the train matrices:
                emg_reduced = emg(setdiff(1:length(chords),ch),:);
                emg_reduced_proj = emg_proj(setdiff(1:length(chords),ch),:);
                emg_reduced_crippled01 = emg_crippled01(setdiff(1:length(chords),ch),:);
                emg_reduced_crippled02 = emg_crippled02(setdiff(1:length(chords),ch),:);
                force_reduced = force(setdiff(1:length(chords),ch),:);
            
                % EMG to force linear mixture matrix:
                W = (emg_reduced' * emg_reduced)^-1 * emg_reduced' * force_reduced;
                W_proj = (emg_reduced_proj' * emg_reduced_proj)^-1 * emg_reduced_proj' * force_reduced;
                W_crippled01 = (emg_reduced_crippled01' * emg_reduced_crippled01)^-1 * emg_reduced_crippled01' * force_reduced;
                W_crippled02 = (emg_reduced_crippled02' * emg_reduced_crippled02)^-1 * emg_reduced_crippled02' * force_reduced;

                % add predictions to the corresponding matrices:
                force_pred = emg(ch,:)*W;
                force_pred_proj = emg_proj(ch,:)*W_proj;
                force_pred_crippled01 = emg_crippled01(ch,:)*W_crippled01;
                force_pred_crippled02 = emg_crippled02(ch,:)*W_crippled02;

                F_crossval(ch,:) = force_pred;
                F_crossval_proj(ch,:) = force_pred_proj;
                F_crossval_crippled01(ch,:) = force_pred_crippled01;
                F_crossval_crippled02(ch,:) = force_pred_crippled02;
            end
            tmp.R2_m2f_crossval = mean(1 - sum((F_crossval-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            tmp.R_m2f_crossval = corr2(F_crossval,force_pattern.pattern{i});  

            tmp.R2_m2f_crossval_proj = mean(1 - sum((F_crossval_proj-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            tmp.R_m2f_crossval_proj = corr2(F_crossval_proj,force_pattern.pattern{i});  

            tmp.R2_m2f_crossval_crippled01 = mean(1 - sum((F_crossval_crippled01-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            tmp.R_m2f_crossval_crippled01 = corr2(F_crossval_crippled01,force_pattern.pattern{i});  

            tmp.R2_m2f_crossval_crippled02 = mean(1 - sum((F_crossval_crippled02-force_pattern.pattern{i}).^2,1) ./ sum(force_pattern.pattern{i}.^2,1));
            tmp.R_m2f_crossval_crippled02 = corr2(F_crossval_crippled02,force_pattern.pattern{i});  

            C = addstruct(C,tmp,'row','force');
        end
        varargout{1} = C;

        % PLOT - Crippled Force Explanation:
        fig = figure('Units','centimeters', 'Position',[15 15 25 25]);
        x = [ones(length(C.R_m2f_crippled01),1) ; 2*ones(length(C.R_m2f_crippled02),1)];
        y = [C.R_m2f_crippled01 ; C.R_m2f_crippled02];
        [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_width); 
        drawline(mean(C.R_m2f),'dir','horz','lim',[0,3.5],'linewidth',4,'color',[0.85 0.85 0.85],'linestyle',':')
        
        ylim([0,1])
        xlim([0,3.5])
        h = gca;
        h.XTickLabels = {num2str(var1,'%.2f'),num2str(var2,'%.2f')};
        h.XAxis.FontSize = my_font.conf_tick_label;
        h.YAxis.FontSize = my_font.conf_tick_label;
        h.LineWidth = 2;
        ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
        fontname("Arial")
        
        % PLOT - var explained by dimensions within chord and natural:
        fig = figure('Units','centimeters', 'Position',[15 15 15 15]);
        drawline(100,'dir','horz','lim',[0,11],'linewidth',4,'color',[0.8 0.8 0.8],'linestyle','-'); hold on;

        x = repmat((1:10)',size(C_PCA.explained_nat,1),1);
        y = C_PCA.explained_nat';
        eig_sum = sum(y,1);
        y = y./eig_sum;
        y = cumsum(y,1)*100;
        y = y(:);
        lineplot(x,y,'markertype','o','markersize',5,'markerfill',colors_red(2,:),'markercolor',colors_red(2,:),'linecolor',colors_red(2,:),'linewidth',2,'errorcolor',colors_red(2,:)); hold on;

        x = repmat((1:10)',size(C_PCA.explained_chord,1),1);        
        y = C_PCA.explained_chord';
        eig_sum = sum(y,1);
        y = y./eig_sum;
        y = cumsum(y,1)*100;
        y = y(:);
        lineplot(x,y,'markersize',5,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',2,'errorcolor',colors_blue(5,:));
        lgd = legend({'','natural','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','chord'});
        xlim([0,11])
        ylim([0,105])
        ylabel('Variance Explained (%)')
        xlabel('Number of Dimensions')

        % PLOT - Projected Chord EMG Variance:
        fig = figure('Units','centimeters', 'Position',[15 15 25 25]);
        x = kron([1:10]',ones(1,length(emg_pattern.pattern)));
        x = x(:);
        y = C.var_chord2chord'./C.total_var_chord' * 100;
        y = y(:);
        lineplot(x,y,'markersize',5,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',2,'errorcolor',colors_blue(5,:)); 
        hold on;
        y = C.var_chord2nat'./C.total_var_chord' * 100;
        y = y(:);
        lineplot(x,y,'markertype','o','markersize',5,'markerfill',colors_red(2,:),'markercolor',colors_red(2,:),'linecolor',colors_red(2,:),'linewidth',2,'errorcolor',colors_red(2,:));
        drawline(100,'dir','horz','lim',[0,11],'linewidth',4,'color',[0.85 0.85 0.85],'linestyle','-')
        
        ylim([0,102])
        xlim([0,11])
        h = gca;
        h.XAxis.FontSize = my_font.conf_tick_label;
        h.YAxis.FontSize = my_font.conf_tick_label;
        h.LineWidth = 2;
        ylabel('%Var of Chord EMG','FontSize',my_font.conf_label)
        xlabel('directions (sorted)','FontSize',my_font.conf_label)
        fontname("Arial")

        % PLOT - Projected Nat EMG Variance:
        fig = figure('Units','centimeters', 'Position',[15 15 25 25]);
        x = kron([1:10]',ones(1,length(emg_pattern.pattern)));
        x = x(:);
        y = C.var_nat2chord'./C.total_var_nat' * 100;
        y = y(:);
        lineplot(x,y,'markersize',5,'markerfill',colors_blue(5,:),'markercolor',colors_blue(5,:),'linecolor',colors_blue(5,:),'linewidth',2,'errorcolor',colors_blue(5,:)); 
        hold on;
        y = C.var_nat2nat'./C.total_var_nat' * 100;
        y = y(:);
        lineplot(x,y,'markertype','o','markersize',5,'markerfill',colors_red(2,:),'markercolor',colors_red(2,:),'linecolor',colors_red(2,:),'linewidth',2,'errorcolor',colors_red(2,:));
        drawline(100,'dir','horz','lim',[0,11],'linewidth',4,'color',[0.85 0.85 0.85],'linestyle','-')
        
        ylim([0,102])
        xlim([0,11])
        h = gca;
        h.XAxis.FontSize = my_font.conf_tick_label;
        h.YAxis.FontSize = my_font.conf_tick_label;
        h.LineWidth = 2;
        ylabel('%Var of Natural EMG','FontSize',my_font.conf_label)
        xlabel('directions (sorted)','FontSize',my_font.conf_label)
        fontname("Arial")
        
    case 'finger_count_explanation'
        C_MD = natChord_analyze('model_testing_all_efc1','measure','MD','model_names',{'n_fing'},'is_plot',0);
        C_RT = natChord_analyze('model_testing_all_efc1','measure','RT','model_names',{'n_fing'},'is_plot',0);
        
        x = [ones(length(C_MD.r),1)];
        y = [C_MD.r];
        fig = figure('Units','centimeters', 'Position',[15 15 12 17.5]);
        [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_line_width,'gapwidth',[0.5,0,0,0],'barwidth',1); 
        ylim([0,1])
        lim_width = 0.5;
        drawline(mean(C_MD.noise_ceil),'dir','horz','lim',[0 2],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
        % xlim([0,6])
        h = gca;
        h.XTickLabels = {'Finger Count Model'};
        h.XAxis.FontSize = my_font.conf_label;
        h.YAxis.FontSize = my_font.conf_tick_label;
        h.YTick = 0:0.2:1;
        h.LineWidth = conf.axis_width;
        ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
        fontname("Arial")
    
        x = [ones(length(C_RT.r),1)];
        y = [C_RT.r];
        fig = figure('Units','centimeters', 'Position',[15 15 12 17.5]);
        [x_coord,~,~] = barplot(x,y,'capwidth',0.1,'linewidth',conf.bar_line_width,'gapwidth',[0.5,0,0,0],'barwidth',1); 
        ylim([0,1])
        lim_width = 0.5;
        drawline(mean(C_RT.noise_ceil),'dir','horz','lim',[0 2],'linewidth',conf.horz_line_width,'color',[0.8 0.8 0.8],'linestyle',':')
        % xlim([0,6])
        h = gca;
        h.YTick = 0:0.5:1;
        h.XTickLabels = {'Finger Count Model'};
        h.XAxis.FontSize = my_font.conf_label;
        h.YAxis.FontSize = my_font.conf_tick_label;
        h.LineWidth = conf.axis_width;
        ylabel('$\mathbf{R}$','interpreter','LaTex','FontSize',my_font.conf_label)
        fontname("Arial")

    case 'PCA_natural'
        sampling_option = 'halves';
        vararginoptions(varargin,{'sampling_option'});
        % loading data:
        [emg_rel,emg_pattern] = natChord_analyze('emg_reliability');
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        
        chords = data.chordID(data.sn==1 & data.sess==1);
        subjects = unique(data.sn);

        % tranforming subject numbers to subject names:
        subject_names = strcat('subj',num2str(subjects,'%02.f'));
        
        % define dataframe:
        C = [];
        % loop on subj:
        for sn = 1:length(subjects)
            % natural EMG file name:
            file_name = fullfile(project_path, 'analysis', ['natChord_' subject_names(sn,:) '_emg_natural_' sampling_option '.mat']);
            
            % loading natural EMG dists:
            emg_dist = load(file_name);
            emg_dist = emg_dist.emg_natural_dist;
            
            sess = unique(data.sess(data.sn == subjects(sn)));
            sess = [1];
            % loop on subject sessions:
            for j = 1:length(sess)
                % scaling factors:
                scales = get_emg_scales(subjects(sn),sess(j));
                
                emg_dist_sess = getrow(emg_dist,emg_dist.sess==sess(j));
                % normalizing the natural EMGs:
                for i = 1:length(emg_dist_sess.dist)
                    emg_dist_sess.dist{i} = (emg_dist_sess.dist{i} ./ scales);
                end
                
                % PCA of half-half natural data:
                half_idx = [1,2];
                for h = 1:2
                    the_half = half_idx(h);
                    natural_emg = emg_dist_sess.dist{the_half};
                    [coeff_nat, score_nat, eigenvalues_nat, ~, explained_nat] = pca(natural_emg);
                    tmp.sn = subjects(sn);
                    tmp.sess = sess(j);
                    tmp.half = h;
                    tmp.explained_nat = explained_nat';
                    tmp.coef_nat = {coeff_nat};
                    tmp.natural_dist = {natural_emg};
                    C = addstruct(C,tmp,'row','force');
                end
            end
        end
        save(fullfile(project_path,'analysis','natural_pca.mat'),"C");
        varargout{1} = C;

    case 'PCA_chord'
        % loading data:
        data = dload(fullfile(project_path, 'analysis', 'natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);
        subjects = unique(data.sn);

        % define dataframe:
        C = [];
        % loop on subj:
        for sn = 1:length(subjects)
            sess = unique(data.sess(data.sn == subjects(sn)));
            % loop on subject sessions:
            for j = 1:length(sess)
                % scaling factors:
                row_data = data.sn==subjects(sn) & data.sess==sess(j);
                emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                scales = get_emg_scales(subjects(sn),sess(j));
                emg = emg ./ scales;
                
                % PCA of half-half EMG data:
                [coeff_chord, score, eigenvalues, ~, explained_var] = pca(emg);
                tmp.sn = subjects(sn);
                tmp.sess = sess(j);
                tmp.explained_var = explained_var';
                tmp.coef_chord = {coeff_chord};
                tmp.emg = {emg};
                C = addstruct(C,tmp,'row','force');
            end
        end
        save(fullfile(project_path,'analysis','chord_pca.mat'),"C");
        varargout{1} = C;

        
    case 'natural_PC_explanation'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        
        subj = unique(C.sn);
        halves  = unique(C.half);
        
        ANA = [];
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            % Ali changed here - only analyze first session
            sess = [1];
            for j = 1:length(sess) 
                % loop on halves:
                for k = 1:length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    coef = C.coef_nat{row_C};
                    row_data = data.sn==subj(i) & data.sess==sess(j);

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    emg = emg - mean(emg,1);
                    M = emg'*emg;
                    emg_var = trace(M);
                    
                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};
                    emg_nat = emg_nat - mean(emg_nat,1);
                    nat = emg_nat' * emg_nat;
                    nat_var = trace(nat);

                    % get projected to PC variance:
                    for i_pc = 1:size(coef,2)
                        var_chord2nat = trace(coef(:,i_pc)' * M * coef(:,i_pc));
                        var_nat2nat = trace(coef(:,i_pc)' * nat * coef(:,i_pc));
                        
                        tmp_df.sn = subj(i);
                        tmp_df.sess = sess(j);
                        tmp_df.half = halves(k);
                        tmp_df.PC = i_pc;
                        tmp_df.nat_explained = var_nat2nat / nat_var * 100;
                        tmp_df.chord_explained = var_chord2nat / emg_var * 100;
                        ANA = addstruct(ANA,tmp_df,'row','force');
                        
                        % if i_pc == 1
                        %     fprintf('\nNat = %.3f , chord=%.3f\n',var_nat2nat / nat_var * 100,var_chord2nat / emg_var * 100)
                        % end
                    end
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_pca.tsv'),ANA);
        varargout{1} = ANA;

    case 'pca_impaired_model'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        
        subj = unique(C.sn);
        halves  = unique(C.half);

        ANA = [];
        % num dimensions to run the emg->force transform:
        dim = 1:size(C.coef_nat{1},2);
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            for j = 1:length(sess)
                % loop on halves:
                for k = 1:length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    row_data = data.sn==subj(i) & data.sess==sess(j);
                    coef = C.coef_nat{row_C};

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    emg = emg - mean(emg,1);
                    M = emg'*emg;
                    emg_var = trace(M);

                    % get the force pattern:
                    % differential forces (force_i = ext_i - flex_i):
                    diff_force = [data.diff_force_f1(row_data),...
                                  data.diff_force_f2(row_data),...
                                  data.diff_force_f3(row_data),...
                                  data.diff_force_f4(row_data),...
                                  data.diff_force_f5(row_data)];
                    % making a 10-D force pattern from the differential:
                    force = [diff_force,diff_force];
                    for col = 1:5
                        force(force(:,col)<0,col) = 0;
                        force(force(:,col+5)>0,col+5) = 0;
                    end
                    force(:,6:10) = force(:,6:10)*-1;
                    force = force - mean(force,1);
                    
                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};
                    emg_nat = emg_nat - mean(emg_nat,1);
                    nat = emg_nat' * emg_nat;
                    nat_var = trace(nat);

                    % running the crippled models:
                    for i_dim = 1:length(dim)
                        % define the dimensions:
                        dim1 = dim(1:i_dim);
                        if (i_dim ~= length(dim))
                            dim2 = dim(i_dim+1:end);
                        else
                            dim2 = [];
                        end

                        % variance calculations:
                        var_nat1 = trace(coef(:,dim1)' * nat * coef(:,dim1)) / nat_var * 100;
                        var_chord1 = trace(coef(:,dim1)' * M * coef(:,dim1)) / emg_var * 100;
                        if ~isempty(dim2)
                            var_nat2 = trace(coef(:,dim2)' * nat * coef(:,dim2)) / nat_var * 100;
                            var_chord2 = trace(coef(:,dim2)' * M * coef(:,dim2)) / emg_var * 100;
                        else
                            var_nat2 = 0;
                            var_chord2 = 0;
                        end

                        % force prediction:
                        emg_crippled01 = emg * coef(:,dim1);
                        emg_crippled01 = emg_crippled01 - mean(emg_crippled01,1);
                        % linear regression:
                        W_crippled01 = (emg_crippled01' * emg_crippled01)^-1 * emg_crippled01' * force;
                        % reconstruct the force:
                        force_pred_crippled01 = emg_crippled01 * W_crippled01;
                        % calculate reconstruction performance:
                        R2_crippled01 = mean(1 - sum((force_pred_crippled01-force).^2,1) ./ sum(force.^2,1));
                        R_crippled01 = corr(force_pred_crippled01(:),force(:));

                        if ~isempty(dim2)
                            emg_crippled02 = emg * coef(:,dim2);
                            emg_crippled02 = emg_crippled02 - mean(emg_crippled02,1);
                            W_crippled02 = (emg_crippled02' * emg_crippled02)^-1 * emg_crippled02' * force;
                            force_pred_crippled02 = emg_crippled02 * W_crippled02;
                            R2_crippled02 = mean(1 - sum((force_pred_crippled02-force).^2,1) ./ sum(force.^2,1));
                            R_crippled02 = corr(force_pred_crippled02(:),force(:));
                        else
                            R2_crippled02 = 0;
                            R_crippled02 = 0;
                        end
                        tmp_df.sn = subj(i);
                        tmp_df.sess = sess(j);
                        tmp_df.half = halves(k);
                        tmp_df.dim1 = length(dim1);
                        tmp_df.dim2 = length(dim2);
                        tmp_df.nat_explained1 = var_nat1;
                        tmp_df.nat_explained2 = var_nat2;
                        tmp_df.chord_explained1 = var_chord1;
                        tmp_df.chord_explained2 = var_chord2;
                        tmp_df.r_force1 = R_crippled01;
                        tmp_df.r_force2 = R_crippled02;
                        tmp_df.r2_force1 = R2_crippled01;
                        tmp_df.r2_force2 = R2_crippled02;
                        ANA = addstruct(ANA,tmp_df,'row','force');
                    end
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_impaired_model.tsv'),ANA);
        varargout{1} = ANA;

    case 'pca_impaired_model_crossvall'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;      
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);

        subj = unique(C.sn);
        halves  = unique(C.half);
        
        ANA = [];
        % num dimensions to run the emg->force transform:
        dim = 1:size(C.coef_nat{1},2);
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            sess = [1]; % Ali changed here.
            for j = 1:length(sess)
                % loop on halves:
                for k = 1:length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    row_data = data.sn==subj(i) & data.sess==sess(j);
                    coef = C.coef_nat{row_C};

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    emg = emg - mean(emg,1);
                    M = emg'*emg;
                    emg_var = trace(M);

                    % get the force pattern:
                    % differential forces (force_i = ext_i - flex_i):
                    diff_force = [data.diff_force_f1(row_data),...
                                  data.diff_force_f2(row_data),...
                                  data.diff_force_f3(row_data),...
                                  data.diff_force_f4(row_data),...
                                  data.diff_force_f5(row_data)];
                    % making a 10-D force pattern from the differential:
                    force = [diff_force,diff_force];
                    for col = 1:5
                        force(force(:,col)<0,col) = 0;
                        force(force(:,col+5)>0,col+5) = 0;
                    end
                    force(:,6:10) = force(:,6:10)*-1;
                    force = force - mean(force,1);
                    
                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};
                    emg_nat = emg_nat - mean(emg_nat,1);
                    nat = emg_nat' * emg_nat;
                    nat_var = trace(nat);

                    % running the crippled models:
                    for i_dim = 1:length(dim)
                        % define the dimensions:
                        dim1 = dim(1:i_dim);
                        dim2 = dim(end-i_dim+1:end);

                        % variance calculations:
                        var_nat1 = trace(coef(:,dim1)' * nat * coef(:,dim1)) / nat_var * 100;
                        var_chord1 = trace(coef(:,dim1)' * M * coef(:,dim1)) / emg_var * 100;
                        var_nat2 = trace(coef(:,dim2)' * nat * coef(:,dim2)) / nat_var * 100;
                        var_chord2 = trace(coef(:,dim2)' * M * coef(:,dim2)) / emg_var * 100;
                        
                        % impaired EMGs:
                        emg_crippled01 = emg * coef(:,dim1);
                        emg_crippled01 = emg_crippled01 - mean(emg_crippled01,1);

                        emg_crippled02 = emg * coef(:,dim2);
                        emg_crippled02 = emg_crippled02 - mean(emg_crippled02,1);
                            
                        % cross-validated force prediction with the
                        % impaired models:
                        F_crossval_crippled01 = zeros(size(force));
                        F_crossval_crippled02 = zeros(size(force));
                        for ch = 1:length(chords)
                            % remove one chord from the train matrices:
                            emg_reduced_crippled01 = emg_crippled01(setdiff(1:length(chords),ch),:);
                            emg_reduced_crippled02 = emg_crippled02(setdiff(1:length(chords),ch),:);
                            force_reduced = force(setdiff(1:length(chords),ch),:);
                        
                            % EMG to force linear mixture matrix:
                            W_crippled01 = (emg_reduced_crippled01' * emg_reduced_crippled01)^-1 * emg_reduced_crippled01' * force_reduced;
                            W_crippled02 = (emg_reduced_crippled02' * emg_reduced_crippled02)^-1 * emg_reduced_crippled02' * force_reduced;
            
                            % add predictions to the corresponding matrices:
                            force_pred_crippled01 = emg_crippled01(ch,:)*W_crippled01;
                            force_pred_crippled02 = emg_crippled02(ch,:)*W_crippled02;
                            
                            F_crossval_crippled01(ch,:) = force_pred_crippled01;
                            F_crossval_crippled02(ch,:) = force_pred_crippled02;
                        end
                        
                        % evaluate the model performance:
                        R2_crippled01 = mean(1 - sum((F_crossval_crippled01-force).^2,1) ./ sum(force.^2,1));
                        R_crippled01 = corr(F_crossval_crippled01(:),force(:));
            
                        R2_crippled02 = mean(1 - sum((F_crossval_crippled02-force).^2,1) ./ sum(force.^2,1));
                        R_crippled02 = corr(F_crossval_crippled02(:),force(:));  

                        tmp_df.sn = subj(i);
                        tmp_df.sess = sess(j);
                        tmp_df.half = halves(k);
                        tmp_df.dim1 = dim1(end);
                        tmp_df.dim2 = dim2(1);
                        tmp_df.nat_explained1 = var_nat1;
                        tmp_df.nat_explained2 = var_nat2;
                        tmp_df.chord_explained1 = var_chord1;
                        tmp_df.chord_explained2 = var_chord2;
                        tmp_df.r_force1 = R_crippled01;
                        tmp_df.r_force2 = R_crippled02;
                        tmp_df.r2_force1 = R2_crippled01;
                        tmp_df.r2_force2 = R2_crippled02;
                        ANA = addstruct(ANA,tmp_df,'row','force');
                    end
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_impaired_model_crossval.tsv'),ANA);
        varargout{1} = ANA;

    case 'chord_pca_impaired_model_crossvall'
        C = load(fullfile(project_path,'analysis','chord_pca.mat')).C;      
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);

        subj = unique(C.sn);
        
        ANA = [];
        % num dimensions to run the emg->force transform:
        dim = 1:size(C.coef_chord{1},2);
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            for j = 1:length(sess)
                % get the coefs:
                row_C = C.sn==subj(i) & C.sess==sess(j);
                row_data = data.sn==subj(i) & data.sess==sess(j);
                coef = C.coef_chord{row_C};

                % get the chord emg patterns:
                emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                       data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                scales = get_emg_scales(subj(i),sess(j));
                emg = emg ./ scales;
                emg = emg - mean(emg,1);
                M = emg'*emg;
                emg_var = trace(M);

                % get the force pattern:
                % differential forces (force_i = ext_i - flex_i):
                diff_force = [data.diff_force_f1(row_data),...
                              data.diff_force_f2(row_data),...
                              data.diff_force_f3(row_data),...
                              data.diff_force_f4(row_data),...
                              data.diff_force_f5(row_data)];
                % making a 10-D force pattern from the differential:
                force = [diff_force,diff_force];
                for col = 1:5
                    force(force(:,col)<0,col) = 0;
                    force(force(:,col+5)>0,col+5) = 0;
                end
                force(:,6:10) = force(:,6:10)*-1;
                force = force - mean(force,1);
                
                % running the crippled models:
                for i_dim = 1:length(dim)
                    % accumulative dimensions from highest order PC to
                    % lowest:
                    dim_forward = dim(1:i_dim);
                    % PC from lowest order to highest order:
                    dim_backward = dim(end-i_dim+1:end);

                    % variance calculations:
                    var_chord1 = trace(coef(:,dim_forward)' * M * coef(:,dim_forward)) / emg_var * 100;
                    var_chord2 = trace(coef(:,dim_backward)' * M * coef(:,dim_backward)) / emg_var * 100;
                    
                    % impaired EMGs:
                    emg_crippled01 = emg * coef(:,dim_forward);
                    emg_crippled01 = emg_crippled01 - mean(emg_crippled01,1);

                    emg_crippled02 = emg * coef(:,dim_backward);
                    emg_crippled02 = emg_crippled02 - mean(emg_crippled02,1);
                        
                    % cross-validated force prediction with the
                    % impaired models:
                    F_crossval_crippled01 = zeros(size(force));
                    F_crossval_crippled02 = zeros(size(force));
                    for ch = 1:length(chords)
                        % remove one chord from the train matrices:
                        emg_reduced_crippled01 = emg_crippled01(setdiff(1:length(chords),ch),:);
                        emg_reduced_crippled02 = emg_crippled02(setdiff(1:length(chords),ch),:);
                        force_reduced = force(setdiff(1:length(chords),ch),:);
                    
                        % EMG to force linear mixture matrix:
                        W_crippled01 = (emg_reduced_crippled01' * emg_reduced_crippled01)^-1 * emg_reduced_crippled01' * force_reduced;
                        W_crippled02 = (emg_reduced_crippled02' * emg_reduced_crippled02)^-1 * emg_reduced_crippled02' * force_reduced;
        
                        % add predictions to the corresponding matrices:
                        force_pred_crippled01 = emg_crippled01(ch,:)*W_crippled01;
                        force_pred_crippled02 = emg_crippled02(ch,:)*W_crippled02;
                        
                        F_crossval_crippled01(ch,:) = force_pred_crippled01;
                        F_crossval_crippled02(ch,:) = force_pred_crippled02;
                    end
                    
                    % evaluate the model performance:
                    R2_crippled01 = mean(1 - sum((F_crossval_crippled01-force).^2,1) ./ sum(force.^2,1));
                    R_crippled01 = corr(F_crossval_crippled01(:),force(:));  
        
                    R2_crippled02 = mean(1 - sum((F_crossval_crippled02-force).^2,1) ./ sum(force.^2,1));
                    R_crippled02 = corr(F_crossval_crippled02(:),force(:));  

                    tmp_df.sn = subj(i);
                    tmp_df.sess = sess(j);
                    tmp_df.dim1 = dim_forward(end);
                    tmp_df.dim2 = dim_backward(1);
                    tmp_df.chord_explained1 = var_chord1;
                    tmp_df.chord_explained2 = var_chord2;
                    tmp_df.r_force1 = R_crippled01;
                    tmp_df.r_force2 = R_crippled02;
                    tmp_df.r2_force1 = R2_crippled01;
                    tmp_df.r2_force2 = R2_crippled02;
                    ANA = addstruct(ANA,tmp_df,'row','force');
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_impaired_chord_model_crossval.tsv'),ANA);
        varargout{1} = ANA;

    case 'single_PC_impaired_crossval'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;    
        C_chord = load(fullfile(project_path,'analysis','chord_pca.mat')).C;   
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);

        subj = unique(C.sn);
        halves  = unique(C.half);
        
        ANA = [];
        % num dimensions to run the emg->force transform:
        dim = 1:size(C.coef_nat{1},2);
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            sess = [1]; % Ali changed here.
            for j = 1:length(sess)
                % loop on halves:
                for k = 1:length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    row_data = data.sn==subj(i) & data.sess==sess(j);
                    % natural PCs:
                    coef = C.coef_nat{row_C};
                    % chord PCs:
                    coef_chord = C_chord.coef_chord{C_chord.sn==subj(i) & C_chord.sess==sess(j)};

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    emg = emg - mean(emg,1);
                    M = emg'*emg;
                    emg_var = trace(M);

                    % get the force pattern:
                    % differential forces (force_i = ext_i - flex_i):
                    diff_force = [data.diff_force_f1(row_data),...
                                  data.diff_force_f2(row_data),...
                                  data.diff_force_f3(row_data),...
                                  data.diff_force_f4(row_data),...
                                  data.diff_force_f5(row_data)];
                    % making a 10-D force pattern from the differential:
                    force = [diff_force,diff_force];
                    for col = 1:5
                        force(force(:,col)<0,col) = 0;
                        force(force(:,col+5)>0,col+5) = 0;
                    end
                    force(:,6:10) = force(:,6:10)*-1;
                    force = force - mean(force,1);
                    
                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};
                    emg_nat = emg_nat - mean(emg_nat,1);
                    nat = emg_nat' * emg_nat;
                    nat_var = trace(nat);

                    % running the crippled models:
                    for i_dim = 1:length(dim)
                        % variance calculations:
                        var_nat_by_nat = trace(coef(:,i_dim)' * nat * coef(:,i_dim)) / nat_var * 100;
                        var_chord_by_nat = trace(coef(:,i_dim)' * M * coef(:,i_dim)) / emg_var * 100;
                        var_nat_by_chord = trace(coef_chord(:,i_dim)' * nat * coef_chord(:,i_dim)) / nat_var * 100;
                        var_chord_by_chord = trace(coef_chord(:,i_dim)' * M * coef_chord(:,i_dim)) / emg_var * 100;
                        
                        % impaired EMGs:
                        emg_crippled = emg * coef(:,i_dim);
                        emg_crippled = emg_crippled - mean(emg_crippled,1);

                        emg_crippled_chord = emg * coef_chord(:,i_dim);
                        emg_crippled_chord = emg_crippled_chord - mean(emg_crippled_chord,1);
                            
                        % cross-validated force prediction with the
                        % impaired models:
                        F_crossval_crippled = zeros(size(force));
                        F_crossval_crippled_chord = zeros(size(force));
                        for ch = 1:length(chords)
                            % remove one chord from the train matrices:
                            emg_reduced_crippled = emg_crippled(setdiff(1:length(chords),ch),:);
                            emg_reduced_crippled_chord = emg_crippled_chord(setdiff(1:length(chords),ch),:);
                            force_reduced = force(setdiff(1:length(chords),ch),:);
                        
                            % EMG to force linear mixture matrix:
                            W_crippled = (emg_reduced_crippled' * emg_reduced_crippled)^-1 * emg_reduced_crippled' * force_reduced;
                            W_crippled_chord = (emg_reduced_crippled_chord' * emg_reduced_crippled_chord)^-1 * emg_reduced_crippled_chord' * force_reduced;
            
                            % add predictions to the corresponding matrices:
                            force_pred_crippled = emg_crippled(ch,:)*W_crippled;    
                            force_pred_crippled_chord = emg_crippled_chord(ch,:)*W_crippled_chord;  

                            F_crossval_crippled(ch,:) = force_pred_crippled;
                            F_crossval_crippled_chord(ch,:) = force_pred_crippled_chord;
                        end
                        
                        % evaluate the model performance:
                        R2_crippled_by_nat = mean(1 - sum((F_crossval_crippled-force).^2,1) ./ sum(force.^2,1));
                        R_crippled_by_nat = corr(F_crossval_crippled(:),force(:));  
                        R2_crippled_by_chord = mean(1 - sum((F_crossval_crippled_chord-force).^2,1) ./ sum(force.^2,1));
                        R_crippled_by_chord = corr(F_crossval_crippled_chord(:),force(:)); 

                        tmp_df.sn = subj(i);
                        tmp_df.sess = sess(j);
                        tmp_df.half = halves(k);
                        tmp_df.dim = i_dim;
                        tmp_df.var_nat_by_nat = var_nat_by_nat;
                        tmp_df.var_chord_by_nat = var_chord_by_nat;
                        tmp_df.var_nat_by_chord = var_nat_by_chord;
                        tmp_df.var_chord_by_chord = var_chord_by_chord;
                        tmp_df.r_force_by_nat = R_crippled_by_nat;
                        tmp_df.r2_force_by_nat = R2_crippled_by_nat;
                        tmp_df.r_force_by_chord = R_crippled_by_chord;
                        tmp_df.r2_force_by_chord = R2_crippled_by_chord;
                        ANA = addstruct(ANA,tmp_df,'row','force');
                    end
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_single_pc_impared.tsv'),ANA);
        varargout{1} = ANA;

    case 'single_PC_force_subject_crossval'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;    
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        chords = data.chordID(data.sn==1 & data.sess==1);

        data_efc1 = dload(fullfile(project_path,'analysis','efc1_chord.tsv'));
        data_efc1 = getrow(data_efc1,ismember(data_efc1.chordID,chords));

        efc1_sn = unique(data_efc1.sn);
        subj = unique(C.sn);
        halves  = unique(C.half);
        
        ANA = [];
        % num dimensions to run the emg->force transform:
        dim = 1:size(C.coef_nat{1},2);
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            for j = 1:length(sess)
                % loop on halves:
                for k = 1:length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    row_data = data.sn==subj(i) & data.sess==sess(j);
                    % natural PCs:
                    coef = C.coef_nat{row_C};

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    emg = emg - mean(emg,1);
                    M = emg'*emg;
                    emg_var = trace(M);

                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};
                    emg_nat = emg_nat - mean(emg_nat,1);
                    nat = emg_nat' * emg_nat;
                    nat_var = trace(nat);
                
                    for i_dim = 1:length(dim)
                        % variance calculations:
                        var_nat_by_nat = trace(coef(:,i_dim)' * nat * coef(:,i_dim)) / nat_var * 100;
                        
                        % impaired EMGs:
                        emg_crippled = emg * coef(:,i_dim);
                        emg_crippled = emg_crippled - mean(emg_crippled,1);
                            
                        R_crippled_by_nat = 0;
                        R2_crippled_by_nat = 0;
                        for i_sn = 1:length(efc1_sn)
                            % get the force pattern:
                            efc_sess = [3,4];
                            % calculate group avg:
                            avg_force = 0;
                            subj_in = setdiff(efc1_sn,efc1_sn(i_sn));
                            for i1 = 1:length(subj_in)
                                % getting the subject session average:
                                tmp_force = [];
                                for j1 = 1:length(efc_sess)
                                    diff_force = [data_efc1.diff_force_f1(data_efc1.sn==subj_in(i1) & data_efc1.sess==efc_sess(j1)),...
                                                  data_efc1.diff_force_f2(data_efc1.sn==subj_in(i1) & data_efc1.sess==efc_sess(j1)),...
                                                  data_efc1.diff_force_f3(data_efc1.sn==subj_in(i1) & data_efc1.sess==efc_sess(j1)),...
                                                  data_efc1.diff_force_f4(data_efc1.sn==subj_in(i1) & data_efc1.sess==efc_sess(j1)),...
                                                  data_efc1.diff_force_f5(data_efc1.sn==subj_in(i1) & data_efc1.sess==efc_sess(j1))];
                                    force = [diff_force,diff_force];
                                    for col = 1:5
                                        force(force(:,col)<0,col) = 0;
                                        force(force(:,col+5)>0,col+5) = 0;
                                    end
                                    force(:,6:10) = force(:,6:10)*-1;
                                    tmp_force = cat(3,tmp_force,force);
                                end
                                tmp_force = sum(tmp_force,3,'omitnan')/length(efc_sess);
                    
                                avg_force = avg_force + tmp_force/length(subj_in);
                            end
                            chordID = data_efc1.chordID(data_efc1.sn==1 & data_efc1.sess==4);
    
                            % making the design matrix:
                            X_train = zeros(length(chords),10);
                            for i1 = 1:length(chords)
                                X_train(i1,:) =  avg_force(chordID==chords(i1),:);
                            end
    
                            % calculate group avg:
                            avg_force = [];
                            for j1 = 1:length(efc_sess)
                                diff_force = [data_efc1.diff_force_f1(data_efc1.sn==efc1_sn(i_sn) & data_efc1.sess==efc_sess(j1)),...
                                              data_efc1.diff_force_f2(data_efc1.sn==efc1_sn(i_sn) & data_efc1.sess==efc_sess(j1)),...
                                              data_efc1.diff_force_f3(data_efc1.sn==efc1_sn(i_sn) & data_efc1.sess==efc_sess(j1)),...
                                              data_efc1.diff_force_f4(data_efc1.sn==efc1_sn(i_sn) & data_efc1.sess==efc_sess(j1)),...
                                              data_efc1.diff_force_f5(data_efc1.sn==efc1_sn(i_sn) & data_efc1.sess==efc_sess(j1))];
                                force = [diff_force,diff_force];
                                for col = 1:5
                                    force(force(:,col)<0,col) = 0;
                                    force(force(:,col+5)>0,col+5) = 0;
                                end
                                force(:,6:10) = force(:,6:10)*-1;
                                avg_force = cat(3,avg_force,force);
                            end
                            avg_force = sum(avg_force,3,'omitnan')/length(efc_sess);
                            % making the design matrix:
                            y = zeros(length(chords),10);
                            for i1 = 1:length(chords)
                                y(i1,:) =  avg_force(chordID==chords(i1),:);
                            end
                            
                            % EMG to force linear mixture matrix:
                            W_crippled = (emg_crippled' * emg_crippled)^-1 * emg_crippled' * X_train;
            
                            % add predictions to the corresponding matrices:
                            y_pred = emg_crippled*W_crippled;
                            
                            % evaluate the model performance:
                            R2_crippled_by_nat = R2_crippled_by_nat + mean(1 - sum((y_pred-y).^2,1) ./ sum(y.^2,1))/length(efc1_sn);
                            R_crippled_by_nat = R_crippled_by_nat + corr(y_pred(:),y(:))/length(efc1_sn);  

                        end
                        tmp_df.sn = subj(i);
                        tmp_df.sess = sess(j);
                        tmp_df.half = halves(k);
                        tmp_df.dim = i_dim;
                        tmp_df.var_nat_by_nat = var_nat_by_nat;
                        tmp_df.r_force_by_nat = R_crippled_by_nat;
                        tmp_df.r2_force_by_nat = R2_crippled_by_nat;
                        ANA = addstruct(ANA,tmp_df,'row','force');
                    end
                end
            end
        end
        dsave(fullfile(project_path,'analysis','natChord_single_PC_force_subject_crossval.tsv'),ANA);
        varargout{1} = ANA;

    case 'avg_chord_EMG_patterns'
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        
        % extension channels:
        e1 = data.emg_hold_avg_e1 ./ data.scale_e1;
        e2 = data.emg_hold_avg_e2 ./ data.scale_e2;
        e3 = data.emg_hold_avg_e3 ./ data.scale_e3;
        e4 = data.emg_hold_avg_e4 ./ data.scale_e4;
        e5 = data.emg_hold_avg_e5 ./ data.scale_e5;
        
        f1 = data.emg_hold_avg_f1 ./ data.scale_f1;
        f2 = data.emg_hold_avg_f2 ./ data.scale_f2;
        f3 = data.emg_hold_avg_f3 ./ data.scale_f3;
        f4 = data.emg_hold_avg_f4 ./ data.scale_f4;
        f5 = data.emg_hold_avg_f5 ./ data.scale_f5;
        
        % claculate the average across subjects and sessions:
        % extension channels:
        [~, ~, e1_avg, chords, ~] = get_sem(e1, ones(size(e1)), ones(size(e1)), data.chordID);
        [~, ~, e2_avg, chords, ~] = get_sem(e2, ones(size(e1)), ones(size(e1)), data.chordID);
        [~, ~, e3_avg, chords, ~] = get_sem(e3, ones(size(e1)), ones(size(e1)), data.chordID);
        [~, ~, e4_avg, chords, ~] = get_sem(e4, ones(size(e1)), ones(size(e1)), data.chordID);
        [~, ~, e5_avg, chords, ~] = get_sem(e5, ones(size(e1)), ones(size(e1)), data.chordID);
        % flexion channels:
        [~, ~, f1_avg, chords, ~] = get_sem(f1, ones(size(f1)), ones(size(f1)), data.chordID);
        [~, ~, f2_avg, chords, ~] = get_sem(f2, ones(size(f1)), ones(size(f1)), data.chordID);
        [~, ~, f3_avg, chords, ~] = get_sem(f3, ones(size(f1)), ones(size(f1)), data.chordID);
        [~, ~, f4_avg, chords, ~] = get_sem(f4, ones(size(f1)), ones(size(f1)), data.chordID);
        [~, ~, f5_avg, chords, ~] = get_sem(f5, ones(size(f1)), ones(size(f1)), data.chordID);

        df = [];
        df.chords = chords;
        df.finger_count = get_num_active_fingers(chords);
        df.e1_avg = e1_avg;
        df.e2_avg = e2_avg;
        df.e3_avg = e3_avg;
        df.e4_avg = e4_avg;
        df.e5_avg = e5_avg;
        df.f1_avg = f1_avg;
        df.f2_avg = f2_avg;
        df.f3_avg = f3_avg;
        df.f4_avg = f4_avg;
        df.f5_avg = f5_avg;
        
        dsave(fullfile(project_path,'analysis','avg_emg_patterns.tsv'), df);

    case 'nmf_test'
        C = load(fullfile(project_path,'analysis','natural_pca.mat')).C;
        data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
        
        subj = unique(C.sn);
        halves  = unique(C.half);

        ANA = [];
        % loop on subj:
        for i = 1:length(subj)
            % loop on sess:
            sess = unique(C.sess(C.sn==subj(i)));
            for j = 1:length(sess)
                % loop on halves:
                for k = 1:1 %length(halves)
                    % get the coefs:
                    row_C = C.sn==subj(i) & C.sess==sess(j) & C.half==halves(k);
                    row_data = data.sn==subj(i) & data.sess==sess(j);

                    % get the chord emg patterns:
                    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
                           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
                    scales = get_emg_scales(subj(i),sess(j));
                    emg = emg ./ scales;
                    
                    % get natural EMG from the other half:
                    emg_nat = C.natural_dist{C.sn==subj(i) & C.sess==sess(j) & C.half==(3-halves(k))};

                    [Wn,Hn,Dn] = nnmf(emg_nat',6);
                    [Wc,Hc,Dc] = nnmf(emg',6);
                end    

                figure;
                for ic = 1:6
                    subplot(6,1,ic)
                    bar(1:10,Wn(:,ic),'k')
                    xlabel('Muscles')
                    sgtitle('nat')
                end
                
                figure;
                for ic = 1:6
                    subplot(6,1,ic)
                    bar(1:10,Wc(:,ic),'k')
                    xlabel('Muscles')
                    sgtitle('chord')
                end
            end
        end

    case 'natural_emg_reliability'
        D = dload(fullfile(baseDir, 'analysis', sprintf('natChord_chord.tsv')));
        subjects = unique(D.sn)';
        
        RDM = {};
        corr_mat = {};
        for sn = subjects
            dist = load(fullfile(baseDir,'analysis',sprintf('natChord_subj%02d_emg_natural_halves.mat',sn)));
            dist = dist.emg_natural_dist;
            
            row = D.sn==sn & D.sess==1;
            emg_scales = [D.scale_e1(row), D.scale_e2(row),D.scale_e3(row),D.scale_e4(row),D.scale_e5(row),...
                          D.scale_f1(row), D.scale_f2(row),D.scale_f3(row),D.scale_f4(row),D.scale_f5(row)];
            emg_scales = emg_scales(1,:);
            % 
            % dist.dist{1} = dist.dist{1} ./ emg_scales;
            % dist.dist{2} = dist.dist{2} ./ emg_scales;
            % 
            % cov1 = dist.dist{1}' * dist.dist{1};
            % cov2 = dist.dist{2}' * dist.dist{2};
            % % cov12 = dist.dist{1}' * dist.dist{2};
            % 
            % figure;
            % imagesc(cov1);
            % 
            % figure;
            % imagesc(cov2);
            fs_force = 500;
            fs_emg = 2148.1481; 
            mov_state = [2,3,4];
            mov_state_names = ['plan','exec','feedback'];
            
            emg = load(fullfile(baseDir,'analysis',sprintf('natChord_subj%02d_emg.mat',sn)));
            emg = emg.EMG_struct;
            mov = load(fullfile(baseDir,'analysis',sprintf('natChord_subj%02d_mov.mat',sn)));
            mov = mov.MOV_struct;
            D_sn = dload(fullfile(baseDir,'analysis',sprintf('natChord_subj%02d_raw.tsv',sn)));
            
            hold_dur = 600; % ms
            chords = [19999, 91999, 99199, 99919, 99991, ...
                      29999, 92999, 99299, 99929, 99992];
            cols_in_emg = [1,2,3,4,5,6,7,8,9,10];
            var_scale = [];
            for j = 1:length(chords)
                chord = chords(j);
                rows = find(D_sn.trialCorr==1 & D_sn.chordID==chord);
                emg_vec = [];
                for i = 1:length(rows)
                    % get mov indices corresponding to hold time
                    idx_end = find(mov{rows(i)}(:,1)==4,1)-1;
                    idx_start = idx_end - round(hold_dur/1000 * fs_force) + 1;
                    
                    % get force:
                    t_f = mov{rows(i)}(idx_start:idx_end,3) - mov{rows(i)}(idx_start,3);
                    force = mov{rows(i)}(idx_start:idx_end,19:23);
                    
                    % find corresponding index in emg data:
                    idx_end_emg = round(mov{rows(i)}(idx_end,3)/1000*fs_emg);
                    idx_start_emg = idx_end_emg - round(hold_dur/1000 * fs_emg) + 1;
                    
                    % select emg data:
                    emg_sig = emg{rows(i)}(idx_start_emg:idx_end_emg,:);
                    t_emg = linspace(0,t_f(end),size(emg_sig,1));
                    
                    % avg_time = mean(emg_sig,1);
                    emg_vec = [emg_vec ; emg_sig(:,cols_in_emg(j))'];
                end
                avg_emg = mean(emg_vec,1);
                err_emg = emg_vec - avg_emg;
                var_scale(j) = var(err_emg(:));
            end
            
            tmp1 = dist.dist{1} ./ sqrt(var_scale);
            tmp1 = dist.dist{1} ./ max(tmp1);
            tmp1 = tmp1 - mean(tmp1,1);
            tmp2 = dist.dist{2} ./ sqrt(var_scale);
            tmp2 = dist.dist{2} ./ max(tmp2);
            tmp2 = tmp2 - mean(tmp2,1);
            
            cov1 = tmp1' * tmp1;
            cov2 = tmp2' * tmp2;
            RDM1 = squareform(pdist(tmp1', 'euclidean'));
            RDM2 = squareform(pdist(tmp2', 'euclidean'));
            
            half1 = [tmp1(1:round(size(tmp1,1)/2),:) ; tmp2(1:round(size(tmp1,1)/2),:)];
            half2 = [tmp1(round(size(tmp1,1)/2)+1:end,:) ; tmp2(round(size(tmp2,1)/2)+1:end,:)];
            corr1 = corr(half1);
            corr2 = corr(half2);
            % corr3 = corr(tmp2(1:round(size(tmp2,1)/2),:));
            % corr4 = corr(tmp2(round(size(tmp2,1)/2)+1:end,:));
            
            % figure;
            % imagesc(sqrt(RDM1));
            % axis square;
            % 
            % figure;
            % imagesc(sqrt(RDM2));
            % axis square;
            % 
            RDM{sn,1} = RDM1;
            RDM{sn,2} = RDM2;

            corr_mat{sn,1} = corr1;
            corr_mat{sn,2} = corr2;
        end
        varargout{1} = RDM;
        varargout{2} = corr_mat;

    otherwise
        error('The analysis you entered does not exist!')
end

