
emg_nat = readtable('data/subj02/emg_natural01.csv');

%%
x = emg_nat;
x(1:3,:) = [];
x = table2array(x);

[row_nan,~] = find(isnan(x));
x(unique(row_nan),:) = [];

emg_data_selected = x(:,2:2:end);
emg_data_selected(:,end) = [];
raw = emg_data_selected;

fs_emg = 2148.1481;         % EMG sampling rate in Hz  
Fstop_lpf = 40;

% EMG filters:
% Define filter parameters
% Design the Butterworth filter
[b, a] = butter(6, Fstop_lpf/(fs_emg/2), 'low'); % 6th order Butterworth filter
% Convert to second-order sections (SOS) format
[sos, g] = tf2sos(b, a);

for j = 1:size(emg_data_selected,2)
    % de-mean rectify EMGs:
    emg_data_selected(:,j) = abs(emg_data_selected(:,j) - mean(emg_data_selected(:,j)));
    
    fprintf("Lowpass Filtering channel %d/%d...\n",j,size(emg_data_selected,2))
    emg_data_selected(:,j) = filtfilt(sos, g, emg_data_selected(:,j));
end


%%
figure;
for i = 1:10
    subplot(5,2,i)
    plot(raw(:,i),'k');
end

figure;
for i = 1:10
    subplot(5,2,i)
    plot(emg_data_selected(:,i),'b')
end




