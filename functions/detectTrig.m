function [riseTimes,riseIdx,fallTimes,fallIdx] = detectTrig(trigSig,timeTrig,emgTime,ampThreshold,numTrials,debugging)
% Ali Ghavampour 2023 - alighavam79@gmail.com
% this function detects the rising edge of the trigger signals 
% from the s626 I/O card. Refer to efc1.cpp project and s626 software docs 
% for details of trigger generation with s626.

% "INPUTS:"
% trigSig: trigger signal
% timeTrig: the time vector for the trig signal
% emgTime: the time vector of EMG signal (It is possible that EMG signal
%          and trigger signal do not have the same sampling freq.
% ampThreshold: the threshold factor for edgde detection
% debugging: 1 or 0 - Turns  on the debugging figures

% This could be used as default values for thresholding the amplitude:
% ampThreshold = 0.4;

% trigger detection
trigSig = trigSig/max(trigSig);     % normalizing the value of trigger.
                                    % IMPORTANT NOTE: 
                                    % you might have to remove the 
                                    % negative sign based on how you
                                    % generate and record your trigger.
                                    % You can simply add a condition to
                                    % automatically determine if the negative
                                    % sign is needed or not. I just wasn't
                                    % much in mood to do it.

invTrig = -trigSig;                 % to find the falling edges.

% hint for another method you could possibly use for trigger detection:
% edgeVec = double(edge(trigSig));   % rising and falling edges
% edgeVec(abs(trigSig)<ampThreshold) = 0;
% edgeIdx = find(edgeVec==1);
% edgeTimes = timeTrig(edgeIdx);

% trigger detection starts here:
diffTrig = diff(trigSig);
diffInvTrig = diff(invTrig);

% thresholding the deviations:
diffTrig(diffTrig < ampThreshold) = 0;
diffInvTrig(diffInvTrig < ampThreshold) = 0;

% Finding the locations:
[pks,locs] = findpeaks(diffTrig);
[pks_inv,locs_inv] = findpeaks(diffInvTrig);

% pks = pks(1:2:end);
% locs = locs(1:2:end);
fprintf("\nNum Trigs Detected = %d , inv %d\n",length(locs),length(locs_inv))
fprintf("Num Trials in Run = %d\n",numTrials)
fprintf("====NumTrial should be equal to NumTrigs====\n\n\n")

% with respect to trigger signal's sampling freq:
riseIdx = locs;
riseTimes = timeTrig(riseIdx);
fallIdx = locs_inv;
fallTimes = timeTrig(fallIdx);

% find the indices with respect to EMG signal's sampling freq:
riseIdx_emg = zeros(size(riseIdx));
fallIdx_emg = zeros(size(fallIdx));
riseTimes_emg = zeros(size(riseTimes));
fallTimes_emg = zeros(size(fallTimes));
for i = 1:length(riseIdx)
    riseIdx_emg(i) = findNearest(emgTime,riseTimes(i));
    fallIdx_emg(i) = findNearest(emgTime,fallTimes(i));
    riseTimes_emg(i) = emgTime(riseIdx_emg(i));
    fallTimes_emg(i) = emgTime(fallIdx_emg(i));
end

% output of the function:
riseTimes = riseTimes_emg;
fallTimes = fallTimes_emg;
riseIdx = riseIdx_emg;
fallIdx = fallIdx_emg;

if (debugging)
    figure;
    hold all
    plot(trigSig,'k','LineWidth',1.5)
    plot(5*diffTrig,'--r','LineWidth',1)
    scatter(locs,pks*5,'red','filled')
    scatter(locs_inv,pks_inv*2,'blue','filled')
    xlabel("time (index)")
    ylabel("Trigger Signal(black), Diff Trigger(red dashed), Detected triggers(red/blue points)")
    ylim([-1.5 1.5])
end

