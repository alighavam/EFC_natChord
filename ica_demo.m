clear

emg_nat = load('analysis/natChord_subj01_emg_natural_whole_sampled.mat');
emg_nat = emg_nat.emg_natural_dist;


sig = emg_nat.dist{1}';

[icasig,A,W] = fastica(sig);


%%

figure;
for i = 1:10
    subplot(5,2,i)
    plot(sig(i,:),'k');
end

figure;
for i = 1:10
    subplot(5,2,i)
    plot(icasig(i,:),'b')
end

