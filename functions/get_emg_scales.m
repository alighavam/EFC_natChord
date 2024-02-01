function scales = get_emg_scales(sn,sess)

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);
project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

data = dload(fullfile(project_path,'analysis','natChord_chord.tsv'));
data = getrow(data, data.sn==sn & data.sess==sess);

scales = [data.scale_e1(1),data.scale_e2(1),data.scale_e3(1),data.scale_e4(1),data.scale_e5(1), ...
                  data.scale_f1(1),data.scale_f2(1),data.scale_f3(1),data.scale_f4(1),data.scale_f5(1)];