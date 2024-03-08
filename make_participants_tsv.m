% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);

participants_tsv = [];

row01 = [];
row01.participant_id = 'subj01';
row01.sex = 'M';
row01.age = 23;
row01.emg_electrode_sess01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row01.emg_electrode_sess02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row01.emg_electrode_nat01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row01.emg_electrode_nat02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';


row02 = [];
row02.participant_id = 'subj02';
row02.sex = 'M';
row02.age = 26;
row02.emg_electrode_sess01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row02.emg_electrode_sess02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row02.emg_electrode_nat01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row02.emg_electrode_nat02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';


row03 = [];
row03.participant_id = 'subj03';
row03.sex = 'M';
row03.age = 26;
row03.emg_electrode_sess01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row03.emg_electrode_sess02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row03.emg_electrode_nat01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row03.emg_electrode_nat02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';

row04 = [];
row04.participant_id = 'subj04';
row04.sex = 'M';
row04.age = 23;
row04.emg_electrode_sess01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row04.emg_electrode_sess02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row04.emg_electrode_nat01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row04.emg_electrode_nat02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';

row05 = [];
row05.participant_id = 'subj05';
row05.sex = 'F';
row05.age = 26;
row05.emg_electrode_sess01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row05.emg_electrode_sess02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row05.emg_electrode_nat01 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';
row05.emg_electrode_nat02 = '7a,6b,5a,5b,6a,10a,8b,7b,8a,10b';


participants_tsv = addstruct(participants_tsv,row01,'row','force');
participants_tsv = addstruct(participants_tsv,row02,'row','force');
participants_tsv = addstruct(participants_tsv,row03,'row','force');
participants_tsv = addstruct(participants_tsv,row04,'row','force');
participants_tsv = addstruct(participants_tsv,row05,'row','force');

dsave(fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord', 'data', 'participants.tsv'),participants_tsv)

