
dat = dload('analysis/natChord_subj01_raw.tsv');
mov = load('analysis/natChord_subj01_mov.mat');
emg = load('analysis/natChord_subj01_emg.mat');
fs_emg = 2148.1481;
fs_force = 500;

finger_count = get_num_active_fingers(dat.chordID);

idx = find(finger_count==1 & dat.trialCorr==1);

t1 = 100;
t2 = 500;