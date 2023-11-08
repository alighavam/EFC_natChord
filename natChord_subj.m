function ANA = natChord_subj(subjName,varargin)

% Handling the input arguments:
smoothing_win_length = 25;  % force smoothing window size in ms
fs_force = 500;             % force signals sampling rate in Hz
fs_emg = 2148.1481;         % EMG sampling rate in Hz   
vararginoptions(varargin,{'smoothing_win_length'});

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);

% set file names:
datFileName = ['data/' subjName '/' 'efc1_' num2str(str2double(subjName(end-1:end))) '.dat'];       % input .dat file
subjFileName = [usr_path '/Desktop/Projects/EFC_natChord/analysis/natChord_' subjName '_raw.tsv'];  % output dat file name (saved in analysis folder)
movFileName = [usr_path '/Desktop/Projects/EFC_natChord/analysis/natChord_' subjName '_mov.mat'];   % output mov file name (saved in analysis folder)
emgFileName = [usr_path '/Desktop/Projects/EFC_natChord/analysis/natChord_' subjName '_emg.mat'];   % output emg file name (saved in analysis folder)
subject_tsv = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord', 'data', subjName);          % tsv file including subject specific parameters

% load .dat file:
D = dload(datFileName);

% container for the dat and mov structs:
ANA = [];
MOV_struct = cell(length(D.BN),1);
EMG_struct = cell(length(D.BN),1);

% EMG filters:
hd = bandpass_filter;
hd_lowpass = lowpass_filter;

oldBlock = -1;
% loop on trials:
for i = 1:length(D.BN)
    % load the mov file for each block:
    if (oldBlock ~= D.BN(i))
        fprintf("Loading the .mov file.\n")
        mov = movload(['data/' subjName '/' 'efc1_' num2str(str2double(subjName(end-1:end))) '_' num2str(D.BN(i),'%02d') '.mov']);
        oldBlock = D.BN(i);
    end
    fprintf('Block: %d , Trial: %d\n',D.BN(i),D.TN(i));
    % trial routine:
    C = efc1_trial(getrow(D,i));

    % adding the routine output to the container:
    ANA = addstruct(ANA,C,'row','force');

    % MOV file: 
    MOV_struct{i} = smoothing(mov{D.TN(i)}, smoothing_win_length, fs);
end

% adding subject name to the struct:
sn = ones(length(D.BN),1) * str2double(subjName(end-1:end));

% remove subNum field:
ANA = rmfield(ANA,'subNum');

% adding subj number to ANA:
ANA.sn = sn;

% saving ANA as a tab delimited file:
dsave(subjFileName,ANA);

% saving mov data as a binary file:
save(movFileName ,'MOV_struct')










