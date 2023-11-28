function sampled_emg = emg_natural_prep(subject_info,emg_data, sess, fs, hd, hd_lpf, wn_type, wn_size, sampling_option, wn_spacing)
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
emg_channels = [{'02'}, emg_channels];

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

    % if channel was a Duo sensor:
    else
        % name of the electrode in table:
        channel_name = ['DuoSensor' emg_channels{i}(1) '_'];
        duo_flag = emg_channels{i}(2);
    end

    % get emg table variable names:
    table_names = emg_data.Properties.VariableNames;
    
    % find the table index corresponding to emg_channels{i}:
    ind = find(~cellfun(@isempty,strfind(table_names,channel_name)));

    % if Avanti:
    if (contains(channel_name,'Avanti'))
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

% filtering the EMG signals:
fprintf("Filtering the raw natural EMG signals:\n\n")
for j = 1:size(emg_data_selected,2)
    fprintf("Bandpass Filtering channel %d/%d...\n",j,size(emg_data_selected,2))
    emg_data_selected(:,j) = abs(filtfilt(hd.Numerator, 1, emg_data_selected(:,j)));
    if (~isempty(hd_lpf))
        fprintf("Lowpass Filtering rectified channel %d/%d...\n",j,size(emg_data_selected,2))
        emg_data_selected(:,j) = filtfilt(hd_lpf.Numerator, 1, emg_data_selected(:,j));
    end
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
            sampled_emg(i,:) = sum(tmp,1)./sum(wn,1);
        end

    case 'whole_sampled'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % selecting the sample intervals:
        intervals = intervals(1:wn_spacing:end,:);

        % sampling the natural EMGs:
        sampled_emg = zeros(size(intervals,1),size(emg_data_selected,2));

        for i = 1:size(intervals,1)
            % windowing the interval:
            tmp = emg_data_selected(intervals(i,1):intervals(i,2),:) .* wn;
            sampled_emg(i,:) = sum(tmp,1)./sum(wn,1);
        end

    case 'whole_thresholded'
        % Creating the sampling intervals:
        intervals = [1:wn_size:size(emg_data_selected,1)-wn_size]';
        intervals = [intervals, [wn_size:wn_size:size(emg_data_selected,1)]'];
        
        % selecting the sample intervals:
        intervals = intervals(1:wn_spacing:end,:);

        % sampling the natural EMGs:
        sampled_emg = zeros(size(intervals,1),size(emg_data_selected,2));

        for i = 1:size(intervals,1)
            % windowing the interval:
            tmp = emg_data_selected(intervals(i,1):intervals(i,2),:) .* wn;
            sampled_emg(i,:) = sum(tmp,1)./sum(wn,1);
        end

        % scaling factors:
        scales = natChord_analyze('get_scale_factor_emg','subject_name',subject_info.participant_id{1});
        
        if (normalize_channels)
            % normalizing the natural EMGs:
            emg_dist{i} = emg_dist{i} ./ scales(:,i)';
        end
        
        % Norm of all samples:
        samples_norm = vecnorm(sampled_emg');

        % avg of the norms:
        avg_norm = mean(samples_norm);

        % sampling the sampled EMGs based on mean norm threshold:
        % indices for each electrode that are more than their avg:
        ind = samples_norm >= avg_norm;

        % sub sampling:
        sampled_emg = sampled_emg(ind,:);

    case 'peaks'

    otherwise
        error('emg_natural_prep: no option %s',sampling_option)
end














