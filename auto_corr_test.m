x = load('analysis/natChord_subj01_emg_natural_whole.mat');

% Assume 'signal' is your input signal, and 'Fs' is the original sampling frequency.
Fs = 1000/20;
signal = x.emg_natural_dist.dist{1};
signal = signal(1:30000,:);

% Assuming signal is a T by C matrix
T = size(signal, 1); % Number of time points
C = size(signal, 2); % Number of channels

% Center the signal by subtracting the mean
signal_centered = signal - mean(signal, 1);

% Compute the multi-variate autocorrelation
max_lag = T - 1; % Set maximum lag
autocorr_values = zeros(max_lag+1, 1); % Initialize the autocorrelation vector

for lag = 0:max_lag
    % Take the product of the signal with its lagged version
    lagged_signal = [zeros(lag, C); signal_centered(1:T-lag, :)];
    
    % Compute the sum of the product across all channels and time points
    autocorr_values(lag + 1) = sum(sum(signal_centered .* lagged_signal)) / (T - lag);
end

% Normalize the autocorrelation values
autocorr_values = autocorr_values / autocorr_values(1);

% Plot the autocorrelation function as a continuous line
figure;
plot(0:max_lag, autocorr_values, '-'); % Using '-' for a solid line
xlabel('Lag');
ylabel('Autocorrelation');
title('Autocorrelation vs Lag for Multivariate Signal');
grid on;
