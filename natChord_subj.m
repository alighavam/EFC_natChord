function ANA = natChord_subj(subjName,varargin)

% Handling the input arguments:
smoothing_win_length = 25;  % force smoothing window size in ms
fs_force = 500;             % force signals sampling rate in Hz
fs_emg = 2148.1481;         % EMG sampling rate in Hz  
Fstop_lpf = 40;
natural_window_size = 100;      % window size to sample natural EMG
sampling_option = 'whole_sampled';      % sampling option to select windows from natural EMGs.
natural_window_type = 'Rect';   % sampling window type for natural EMGs.
wn_spacing = 4;                 % sampling spacing for the 'whole_sampled' option.
vararginoptions(varargin,{'smoothing_win_length','lpf','Fpass_lpf','Fstop_lpf', ...
                          'sampling_option','natural_window_size','natural_window_type','wn_spacing'});

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);

project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

% set file names:
datFileName = fullfile(project_path, 'data', subjName, ['efc1_', num2str(str2double(subjName(end-1:end))), '.dat']); % input .dat file
subjFileName = fullfile(project_path, 'analysis', ['natChord_' subjName '_raw.tsv']);% output dat file name (saved in analysis folder)
movFileName = fullfile(project_path, 'analysis', ['natChord_' subjName '_mov.mat']); % output mov file name (saved in analysis folder)
emgFileName = fullfile(project_path, 'analysis', ['natChord_' subjName '_emg.mat']); % output emg file name (saved in analysis folder)
participants_tsv = fullfile(project_path, 'data', 'participants.tsv');     % tsv file including subject specific parameters

% load participants.tsv:
subj_info = dload(participants_tsv);

% load .dat file:
D = dload(datFileName);
sess = (D.BN<=10) + 2*(D.BN>=11);
sess = unique(sess);

% container for the dat and mov structs:
ANA = [];
MOV_struct = cell(length(D.BN),1);
EMG_struct = cell(length(D.BN),1);

% EMG filters:
% Define filter parameters
% Design the Butterworth filter
[b, a] = butter(6, Fstop_lpf/(fs_emg/2), 'low'); % 6th order Butterworth filter
% Convert to second-order sections (SOS) format
[sos, g] = tf2sos(b, a);

oldBlock = -1;
% loop on trials:
for i = 1:length(D.BN)
    % load the mov file of the block:
    if (oldBlock ~= D.BN(i))
        fprintf("Loading the .mov file. Block %d\n",D.BN(i))
        mov = movload(fullfile(project_path, 'data', subjName, ['efc1_' num2str(str2double(subjName(end-1:end))) '_' num2str(D.BN(i),'%02d') '.mov']));
        
        % load the emg file of the block (this is chord EMG not natural,
        % we'll deal with the natural EMGs later in the code):
        fprintf("Loading emg file %d.\n",D.BN(i))
        emg_data = readtable(fullfile(project_path,'data', subjName, ['emg_run' num2str(D.BN(i),'%02d') '.csv']));
        
        % call emg_chord_prep function:
        [emg_block,baseline_emg,hold_avg_EMG] = emg_chord_prep(emg_data, fs_emg, getrow(D,D.BN == D.BN(i)), mov, ...
                                                getrow(subj_info,find(strcmp(subj_info.participant_id,subjName))), ...
                                                sos,g);
        
        oldBlock = D.BN(i);
    end
    fprintf('Block: %d , Trial: %d\n',D.BN(i),D.TN(i));
    
    % trial routine:
    C = natChord_trial(getrow(D,i),baseline_emg(D.TN(i),:),hold_avg_EMG(D.TN(i),:));

    % adding the trial routine output to the container:
    ANA = addstruct(ANA,C,'row','force');

    % MOV file: 
    MOV_struct{i} = smoothing(mov{D.TN(i)}, smoothing_win_length, fs_emg);
    
    % EMG file:
    EMG_struct{i} = emg_block{D.TN(i)};
end

% adding subject name to the struct:
sn = ones(length(D.BN),1) * str2double(subjName(end-1:end));

% remove subNum field:
ANA = rmfield(ANA,'subNum');

% adding subj number to ANA:
ANA.sn = sn;

% saving ANA as a tab delimited file:
dsave(subjFileName,ANA);

% saving mov and EMG data as binary files:
save(movFileName, 'MOV_struct')
save(emgFileName, 'EMG_struct')

% Preprocessing and dealing with the natural EMGs:
fprintf("Processing natural EMG data...\n\n")

% give session number to the function to know how to create the natural 
% dist:
sess_cell = cellfun(@(x) ['sess', sprintf('%02d', x)], num2cell(sess), 'UniformOutput', false);
make_natural_emg(subjName,sess_cell,fs_emg,hd,hd_lpf,natural_window_type,natural_window_size,sampling_option,wn_spacing);




















