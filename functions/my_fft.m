function [f, P1] = my_fft(x, fs)
% my_fft calculates the single sided spectrum of the input signals.
%
%   Inputs:
%       x:  is a vector or matrix. In case of a matrix, the functions
%       treats each column of x as a vector and calculates the spectrum for
%       each column.
%       
%       fs: sampling frequency of the signals in Hz.
%
%   Outputs:
%       P1:  one sided spectrum of the input signals. If x is a matrix, y
%       will be a matrix of the equal size where each column is the
%       spectrum of the corresponding column in x. The output length equals
%       to the next power of 2 of the original signal length.
%       
%       f:  the frequency domain containing frequencies from 0 to fs/2. The
%       resolution of the frequency domain increases with the length of the
%       signal.

% length of the signal:
if (isvector(x))
    L = length(x);
else
    L = size(x,1);
end

% fft algorithm performs significantly faster if the length of signal is a
% factor of 2. Calculating the next power of 2 from the length of the
% original signal:
n = 2^nextpow2(L);

% fft of x:
Y = fft(x,n);

% two sided spectrum of x:
P2 = abs(Y/n);

% one sided spectrum of x:
P1 = P2(1:n/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);

% defining the frequency domain:
f = fs/n*(0:(n/2));

% plotting:
figure;
plot(f,P1(:,1),"LineWidth",2,'Color','k') 
title("Single-Sided Amplitude Spectrum")
xlabel("f (Hz)")
ylabel("|P1(f)|")



