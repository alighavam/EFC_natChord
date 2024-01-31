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


% row02 = [];
% row02.participant_id = 'subj02';
% row02.sex = 'M';
% row02.age = 25;
% row02.emg_electrode_sess01 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row02.emg_electrode_sess02 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row02.emg_electrode_nat01 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row02.emg_electrode_nat02 = '04,03,13,01,09,5a,10,6b,6a,5b';
% 
% 
% row03 = [];
% row03.participant_id = 'subj03';
% row03.sex = 'F';
% row03.age = 26;
% row03.emg_electrode_sess01 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row03.emg_electrode_sess02 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row03.emg_electrode_nat01 = '04,03,13,01,09,5a,10,6b,6a,5b';
% row03.emg_electrode_nat02 = '04,03,13,01,09,5a,10,6b,6a,5b';


participants_tsv = addstruct(participants_tsv,row01,'row','force');
% participants_tsv = addstruct(participants_tsv,row02,'row','force');
% participants_tsv = addstruct(participants_tsv,row03,'row','force');

dsave(fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord', 'data', 'participants.tsv'),participants_tsv)

