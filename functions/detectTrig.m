function [riseTimes,riseIdx,fallTimes,fallIdx] = detectTrig(trigSig,timeTrig,ampThreshold,numTrials,debugging)
% Ali Ghavampour 2023 - alighavam79@gmail.com
% this function detects the rising edge of the trigger signals 
% from the finger box s626. Refer to efc1.cpp and s626 software docs 
% for details of trigger generation with s626 
% "Inputs help:"
% trig: trigger signal
% tTrig: the time vector for the trig signal
% thresholdFactor: the threshold factor for edgde detection
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

invTrig = -trigSig;             % to find the falling edges.

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

riseIdx = locs;
riseTimes = timeTrig(riseIdx);

fallIdx = locs_inv;
fallTimes = timeTrig(fallIdx);

if (debugging)
    figure;
    hold all
    plot(trigSig,'k','LineWidth',1.5)
    plot(diffTrig,'--r','LineWidth',1)
    scatter(locs,pks,'red','filled')
    scatter(locs_inv,-pks_inv,'blue','filled')
    xlabel("time (index)")
    ylabel("Trigger Signal(black), Diff Trigger(red dashed), Detected triggers(red/blue points)")
    ylim([-1.5 1.5])

end

