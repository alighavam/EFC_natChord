function dist = emg_natural_prep(subject_info,emg_data, sess, fs, sos, g, wn_type, wn_size, sampling_option, wn_spacing)
% Description:
%   wn_size is in ms. 
%   wn_type is 'Gaussian' or 'Rect'
%
%   sampling options: 'whole':  
%                           samples all the non-overlapping windows of 
%                           size window_size from the natural EMG.
%
%                     'whole_sampled': 
%                           First do the 'whole' sampling. Then select
%                           every n samples
%
%                     'whole_thresholded': 
%                           First do the 'whole' sampling. Then, make a
%                           distribution of samples for each channel.
%                           Select some samples for each channel that
%                           follow a selection criterion. For now, we
%                           choose all the windows that are more than the
%                           avg of all windows of a channel. Also, the
%                           windows that are selected are non-overlapping.
%                           So, if something is important for more than one
%                           channel it is not selected twice.
%                           
%                     'peaks':  
%                           samples non-overlapping windows of size
%                           window_size around the local peaks of the
%                           natural EMG.


% The first 3 rows of the loaded emg files is nan:
emg_data(1:3,:) = [];

% if first session:
if (strcmp(sess,'sess01'))
    % getting the used emg channels from subject_info structure:
    emg_channels = strsplit(subject_info.emg_electrode_nat01{1},',');
else
    % getting the used emg channels from subject_info structure:
    emg_channels = strsplit(subject_info.emg_electrode_nat02{1},',');
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

% removing trigger signal from EMG data:
emg_data_selected(:,1:2) = [];

% getting rid of the time vectors from the EMG data:
emg_data_selected = 1000*emg_data_selected(:,2:2:end);

% removing extra NaN elements from the end of the EMG signals (NaNs are
% there bc the sampling freq of trigger channel is not the same as EMG
% channel (BTW THANK YOU "DELSYS" FOR THIS DISCREPANCY AND CHARGING 4 GRANDS 
% PER ELECTRODE :~o ):
[row_nan,~] = find(isnan(emg_data_selected));
emg_data_selected(unique(row_nan),:) = [];

% filtering the EMG signals:
fprintf("Filtering the raw natural EMG signals:\n\n")
for j = 1:size(emg_data_selected,2)
    % de-mean rectify EMGs:
    emg_data_selected(:,j) = abs(emg_data_selected(:,j) - mean(emg_data_selected(:,j)));
    
    fprintf("Lowpass Filtering channel %d/%d...\n",j,size(emg_data_selected,2))
    emg_data_selected(:,j) = filtfilt(sos, g, emg_data_selected(:,j));
end

% converting window size from ms to index:
wn_size = 2*floor(wn_size/2 / 1000 * fs)+1;

% creating the window:
if (strcmp(wn_type,'Rect'))
    wn = ones(wn_size,size(emg_data_selected,2));
elseif (strcmp(wn_type,'Gaussian'))
    wn = gausswin(wn_size,3);
    wn = repmat(wn,1,size(emg_data_selected,2));
else
    error('emg_natural_prep: wn_type %s does not exist.',wn_type)
end


% dataframe to hold dist:
dist = [];
switch sampling_option
    case 'whole'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];

        % sampling the natural EMGs:
        sampled_emg = zeros(size(intervals,1),size(emg_data_selected,2));

        for i = 1:size(intervals,1)
            % windowing the interval:
            tmp = emg_data_selected(intervals(i,1):intervals(i,2),:) .* wn;
            % RMS of window:
            sampled_emg(i,:) = sqrt(sum(tmp.^2,1)./sum(wn,1));
        end
        dist.partition(1,1) = 1;
        dist.dist{1,1} = sampled_emg;

    case 'whole_sampled'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % loop on partitions:
        for i = 1:wn_spacing
            % selecting the sample intervals:
            intervals_tmp = intervals(i:wn_spacing:end,:);

            % sampling the natural EMGs:
            sampled_emg = zeros(size(intervals_tmp,1),size(emg_data_selected,2));
    
            for j = 1:size(intervals_tmp,1)
                % windowing the interval:
                tmp = emg_data_selected(intervals_tmp(j,1):intervals_tmp(j,2),:) .* wn;
                sampled_emg(j,:) = sqrt(sum(tmp.^2,1)./sum(wn,1));
            end
            
            dist.sess(i,1) = str2double(sess(end-1:end));
            dist.partition(i,1) = i;
            dist.dist{i,1} = sampled_emg;
        end

    case 'whole_thresholded'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % loop on partitions:
        for i = 1:wn_spacing
            % selecting the sample intervals:
            intervals_tmp = intervals(i:wn_spacing:end,:);

            % sampling the natural EMGs:
            sampled_emg = zeros(size(intervals_tmp,1),size(emg_data_selected,2));
    
            for j = 1:size(intervals_tmp,1)
                % windowing the interval:
                tmp = emg_data_selected(intervals_tmp(j,1):intervals_tmp(j,2),:) .* wn;
                sampled_emg(j,:) = sqrt(sum(tmp.^2,1)./sum(wn,1));
            end

            % thresholding based on norm of channels after scaling:
            tmp_sample = sampled_emg;

            % scaling factors:
            subject_name = subject_info.participant_id{1};
            scales = get_emg_scales(str2double(subject_name(end-1:end)),str2double(sess(end-1:end)));
            
            % normalizing the natural EMGs:
            tmp_sample = tmp_sample ./ scales;
            
            % Norm of all samples:
            samples_norm = vecnorm(tmp_sample');
    
            % avg of the norms:
            avg_norm = mean(samples_norm);
    
            % sampling the sampled EMGs based on mean norm threshold:
            % indices for each electrode that are more than their avg:
            ind = samples_norm >= 2*avg_norm;
    
            % sub sampling:
            sampled_emg = sampled_emg(ind,:);

            dist.sess(i,1) = str2double(sess(end-1:end));
            dist.partition(i,1) = i;
            dist.dist{i,1} = sampled_emg;
        end

    case 'halves'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % sampling the natural EMGs:
        sampled_emg = zeros(size(intervals,1),size(emg_data_selected,2));

        for i = 1:size(intervals,1)
            % windowing the interval:
            tmp = emg_data_selected(intervals(i,1):intervals(i,2),:) .* wn;
            % RMS of window:
            sampled_emg(i,:) = sqrt(sum(tmp.^2,1)./sum(wn,1));
        end
        dist.sess(1,1) = str2double(sess(end-1:end));
        dist.partition(1,1) = 1;
        dist.sess(2,1) = str2double(sess(end-1:end));
        dist.partition(2,1) = 2;
        
        quarters = round(linspace(1,size(sampled_emg,1),5));
        q1_idx = [quarters(1):quarters(2)];
        q2_idx = [quarters(2)+1:quarters(3)];
        q3_idx = [quarters(3)+1:quarters(4)];
        q4_idx = [quarters(4)+1:quarters(5)];

        dist.dist{1,1} = sampled_emg([q1_idx,q3_idx],:);
        dist.dist{2,1} = sampled_emg([q2_idx,q4_idx],:);

    case 'halves_thresholded'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % sampling the natural EMGs:
        sampled_emg = zeros(size(intervals,1),size(emg_data_selected,2));

        for i = 1:size(intervals,1)
            % windowing the interval:
            tmp = emg_data_selected(intervals(i,1):intervals(i,2),:) .* wn;
            % RMS of window:
            sampled_emg(i,:) = sqrt(sum(tmp.^2,1)./sum(wn,1));
        end
        
        % thresholding based on norm of channels after scaling:
        tmp_sample = sampled_emg;
        
        % scaling factors:
        subject_name = subject_info.participant_id{1};
        scales = get_emg_scales(str2double(subject_name(end-1:end)),str2double(sess(end-1:end)));
        
        % normalizing the natural EMGs:
        tmp_sample = tmp_sample ./ scales;
            
        % Norm of all samples:
        samples_norm = vecnorm(tmp_sample');

        % avg of the norms:
        avg_norm = mean(samples_norm);
    
        % sampling the sampled EMGs based on mean norm threshold:
        % indices for each electrode that are more than their avg:
        ind = samples_norm >= 2*avg_norm;
    
        % sub sampling:
        sampled_emg = sampled_emg(ind,:);

        dist.sess(1,1) = str2double(sess(end-1:end));
        dist.partition(1,1) = 1;
        dist.sess(2,1) = str2double(sess(end-1:end));
        dist.partition(2,1) = 2;
        
        quarters = round(linspace(1,size(sampled_emg,1),5));
        q1_idx = [quarters(1):quarters(2)];
        q2_idx = [quarters(2)+1:quarters(3)];
        q3_idx = [quarters(3)+1:quarters(4)];
        q4_idx = [quarters(4)+1:quarters(5)];

        dist.dist{1,1} = sampled_emg([q1_idx,q3_idx],:);
        dist.dist{2,1} = sampled_emg([q2_idx,q4_idx],:);
    otherwise
        error('emg_natural_prep: no option %s',sampling_option)
end














