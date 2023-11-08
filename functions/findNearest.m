function [nearestIdx] = findNearest(vec,num)

[~,idx] = min(abs(vec-num));
minVal = vec(idx);
nearestIdx = find(vec == minVal,1);