function [ snr ] = SNR( fftV )
%SNR Summary of this function goes here
%   Detailed explanation goes here
[v,vloc] = findpeaks(fftV,'SortStr','descend');
S = v(1).^2;

N = sum(v.^2) - S;
snr = sqrt(S/N);

end

