function sampled_emg = make_natural_emg(subjName,sess,fs_emg,hd,hd_lpf, wn_type, wn_size, sampling_option, wn_spacing)

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);

project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

% output natural emg file name (saved in analysis folder):
nat_emgFileName = fullfile(project_path, 'analysis', ['natChord_' subjName '_emg_natural_' sampling_option '.mat']);

% tsv file including subject specific parameters:
participants_tsv = fullfile(project_path, 'data', 'participants.tsv');     

% load participants.tsv:
subj_info = dload(participants_tsv);

% looping through raw natural emg files:
emg_natural_dist = [];
for i = 1:length(sess)
    fprintf('loading raw natural EMG file %d\n\n',i)
    
    % loading raw natural EMGs:
    emg_data = readtable(fullfile(project_path,'data', subjName, ['emg_natural' num2str(i,'%02d') '.csv']));

    % call emg_natural_prep function:
    sampled_emg = emg_natural_prep(getrow(subj_info,find(strcmp(subj_info.participant_id,subjName))), emg_data, sess{i}, fs_emg, hd, hd_lpf, ...
        wn_type, wn_size, sampling_option, wn_spacing);
    
    emg_natural_dist = addstruct(emg_natural_dist,sampled_emg,'row','force');
end

save(nat_emgFileName,'emg_natural_dist')




