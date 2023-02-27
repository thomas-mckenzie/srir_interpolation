function [hmin] = designMinPhase(H, varargin)
% creates a minimum phase FIR Filter to a given 
% magnitude frequency response using cepstral windowing.
% cepstral window can be adjusted
% inputs:   H ... Magnitude Frequency Responses for positive 
%                 Frequencies dim(H) = (nfft/2+1) x D
%           cepsWindowFraction ... default is 2, increase makes filter
%           even shorter
% nmk20

nfft = (size(H, 1) - 1) * 2; 

cwf = 2; 

if nargin > 1
    cwf = varargin{1}; 
end

H = abs(H); 
H = [H; conj(flipud(H(2:end-1, :)))];
ceps = ifft(log(H));
w = [1; 2*ones(nfft/cwf-1,1); 1; zeros(nfft/cwf-1,1)];
Hmin = exp(fft(ceps.*w));
hmin = ifft(Hmin);

end

