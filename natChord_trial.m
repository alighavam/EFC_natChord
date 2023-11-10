function C = natChord_trial(row,baseline_emg,hold_avg_EMG)

% container for the trial data:
C = [];

C = addstruct(C,row,'row','force');

% adding EMG baseline avg and hold avg to the subject struct:
C.emg_baseline_e1 = baseline_emg(1);
C.emg_baseline_e2 = baseline_emg(2);
C.emg_baseline_e3 = baseline_emg(3);
C.emg_baseline_e4 = baseline_emg(4);
C.emg_baseline_e5 = baseline_emg(5);
C.emg_baseline_f1 = baseline_emg(6);
C.emg_baseline_f2 = baseline_emg(7);
C.emg_baseline_f3 = baseline_emg(8);
C.emg_baseline_f4 = baseline_emg(9);
C.emg_baseline_f5 = baseline_emg(10);


C.emg_hold_avg_e1 = hold_avg_EMG(1);
C.emg_hold_avg_e2 = hold_avg_EMG(2);
C.emg_hold_avg_e3 = hold_avg_EMG(3);
C.emg_hold_avg_e4 = hold_avg_EMG(4);
C.emg_hold_avg_e5 = hold_avg_EMG(5);
C.emg_hold_avg_f1 = hold_avg_EMG(6);
C.emg_hold_avg_f2 = hold_avg_EMG(7);
C.emg_hold_avg_f3 = hold_avg_EMG(8);
C.emg_hold_avg_f4 = hold_avg_EMG(9);
C.emg_hold_avg_f5 = hold_avg_EMG(10);