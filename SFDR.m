function [val] = SFDR(fftR)
% locate max indices
% assumes signal > harmonics
% ingores DC
[v,loc] = findpeaks(fftR,'SortStr','descend');
val = v(1)/v(2);
end