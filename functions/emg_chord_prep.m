function [emg,baseline_emg,hold_avg_EMG] = emg_chord_prep(emg_data, fs, D_block, mov, subject_info, hd, hd_lpf)

% The first 3 rows of the loaded emg files is nan:
emg_data(1:3,:) = [];

% if first session:
if (D_block.BN(1) <= 10)
    % getting the used emg channels from subject_info structure:
    emg_channels = strsplit(subject_info.emg_electrode_sess01{1},',');
else
    % getting the used emg channels from subject_info structure:
    emg_channels = strsplit(subject_info.emg_electrode_sess02{1},',');
end
% sanity check:
if length(emg_channels) ~= 10
    error('EMG channels read from participants.tsv does not have 10 channels! participantID: %s', ...
        subject_info.participant_id{1})
end

% adding trigger channel to emg_channels:
emg_channels = [{'AnalogInputAdapter'}, emg_channels];

% building electrode names based on channels loaded from participant.tsv:
% The container to store the data we select from the emg_table:
emg_data_selected = [];
for i = 1:length(emg_channels)
    % finding channels number:
    channel_num = str2double(emg_channels{i});

    % if channel was an Avanti sensor:
    if (~isnan(channel_num))
        % name of the electrode in table:
        channel_name = ['AvantiSensor' num2str(channel_num) '_'];
    
    % if channel was analog trigger adapter:
    elseif (contains(emg_channels{i},'AnalogInput'))
        channel_name = 'AnalogInputAdapter';

    % if channel was a Duo sensor:
    else
        % name of the electrode in table:
        if length(emg_channels{i})==2
            channel_name = ['DuoSensor' emg_channels{i}(1) '_'];
            duo_flag = emg_channels{i}(2);
        elseif length(emg_channels{i})==3
            channel_name = ['DuoSensor' emg_channels{i}(1:2) '_'];
            duo_flag = emg_channels{i}(3);
        end
    end

    % get emg table variable names:
    table_names = emg_data.Properties.VariableNames;
    
    % find the table index corresponding to emg_channels{i}:
    ind = find(contains(table_names,channel_name));
    
    % if analog tigger:
    if (contains(channel_name,'AnalogInputAdapter'))
        emg_data_selected = [emg_data_selected, table2array(emg_data(:,ind:ind+1))];

    % if Avanti:
    elseif (contains(channel_name,'Avanti'))
        emg_data_selected = [emg_data_selected, table2array(emg_data(:,ind:ind+1))];

    % if Duo:
    else
        if (duo_flag == 'a')
            emg_data_selected = [emg_data_selected, table2array(emg_data(:,ind:ind+1))];
        else
            emg_data_selected = [emg_data_selected, table2array(emg_data(:,ind+2:ind+3))];
        end
    end
    
end

% Extracting triggers of the emg:
t_trig = emg_data_selected(:,1);
t_emg = emg_data_selected(:,3);
trig = emg_data_selected(:,2);

% removing trigger signal from EMG data:
emg_data_selected(:,1:2) = [];

% detecting trial start times from the EMGs:
[~,riseIdx,~,fallIdx] = detectTrig(trig,t_trig,t_emg,0.06,length(D_block.BN),1);

% if number of triggers did not make sense, throw and error:
if (length(riseIdx) ~= length(fallIdx) || length(riseIdx) ~= length(D_block.BN))
    error('Trigger detection went wrong! Block Number = %d',D_block.BN(1))
end

% getting rid of the time vectors from the EMG data:
emg_data_selected = emg_data_selected(:,2:2:end);

% removing extra NaN elements from the end of the EMG signals (NaNs are
% there bc the sampling freq of trigger channel is not the same as EMG
% channel (BTW THANK YOU "DELSYS" FOR THIS DISCREPANCY AND CHARGING 4 GRANDS 
% PER ELECTRODE :~o ):
[row_nan,~] = find(isnan(emg_data_selected));
emg_data_selected(unique(row_nan),:) = [];

% filtering the EMG signals:
fprintf("Filtering the raw EMG signals:\n\n")
for j = 1:size(emg_data_selected,2)
    fprintf("Bandpass Filtering channel %d/%d...\n",j,size(emg_data_selected,2))
    emg_data_selected(:,j) = abs(filtfilt(hd.Numerator, 1, emg_data_selected(:,j)));
    if (~isempty(hd_lpf))
        fprintf("Lowpass Filtering rectified channel %d/%d...\n",j,size(emg_data_selected,2))
        emg_data_selected(:,j) = filtfilt(hd_lpf.Numerator, 1, emg_data_selected(:,j));
    end
end

emg = cell(length(riseIdx),1);
baseline_emg = zeros(length(riseIdx),size(emg_data_selected,2));
hold_avg_EMG = zeros(length(riseIdx),size(emg_data_selected,2));
subset_D = D_block;
for j = 1:length(riseIdx)
    % start of a trial:
    idx1 = riseIdx(j);

    % end of a trial:
    idx2 = fallIdx(j);
    
    % extracting EMG of the trial. Should be T(time) by 10(number of electrodes)
    emg_trial = 1000*emg_data_selected(idx1:idx2,:);
    
    % ===========================================================================
    % Should become a function:
    % ===========================================================================
    % the index of mov file that go cue is presented:
    go_cue_index = find(mov{j}(:,1) == 3, 1);
    
    % the time of that index:
    go_cue_time = mov{j}(go_cue_index,3);

    % planning start time:
    plan_onset_time = (go_cue_time - D_block.planTime(1))/1000;
    
    % finding the plan time start index with respect to EMG:
    idx_plan_onset = round(plan_onset_time * fs);

    % If the beginning of the planning time was at t=0, the index
    % should be 1 because MATLAB's first index is 1:
    if (idx_plan_onset == 0)
        idx_plan_onset = 1;
    end

    fprintf('\n\n')

    baseline_avg_tmp = mean(emg_trial(idx_plan_onset:idx_plan_onset+round(D_block.planTime(1)/1000*fs),:), 1);
    % ===========================================================================
    
    % storing the baseline values:
    baseline_emg(j,:) = baseline_avg_tmp;

    % ===========================================================================
    % Should become a function:
    % ===========================================================================
    % if trial is correct:
    if (subset_D.trialCorr(j))
        % the index that subject achieves the chord -> holds for 600ms:
        end_exec_idx = find(mov{j}(:,1) == 4, 1);
        
        % the time of that index:
        end_exec_time = mov{j}(end_exec_idx,3);

        % hold onset time:
        hold_onset_time = (end_exec_time - 600)/1000;
        
        % finding the hold onset index with respect to EMG:
        idx_hold_onset = round(hold_onset_time * fs);
        
        % averaging EMG during hold time:
        hold_avg_tmp = mean(emg_trial(idx_hold_onset:end,:), 1);
    else
        % if trial is incorrect, we don't have a hold time, hence:
        hold_avg_tmp = -1;
    end

    % ===========================================================================

    % storing the avg hold time EMGs:
    hold_avg_EMG(j,:) = hold_avg_tmp;

    % adding trials to emg cell:
    emg{j} = emg_trial;
end
