function [filtLo,filtHi,fcHz] = ambisonic_crossover(ambisonic_order_or_fc,Fs)
%{
Function to produce linear-phase crossover network for dual-band Ambisonic 
decoding.

Thomas McKenzie, University of York, 2019.

References:
See Thomas McKenzie PhD thesis. 
%}

if ambisonic_order_or_fc <= 49 % if the first argument is a value of 49 or smaller, it is assumed the input is the ambisonic order. Otherwise, it is assumed the input is the cut off frequency. 
    R = 0.09; % radius of reproduction area in metres (Neumann KU 100 is stated on the website as having a width of 18cm)
    fcHz = (ambisonic_order_or_fc * 343) / (4*R*(ambisonic_order_or_fc + 1) * sin(pi / (2*ambisonic_order_or_fc+2))); % this is equation (2) from \cite{Bertet2013}
else
    fcHz = ambisonic_order_or_fc;
end
fcNorm = fcHz/(Fs/2); % convert to normalised frequency

XoverOrder = 128; % crossover order
ripple1 = 50; % ripple in dB

% using chebyshev windows
filtLo = fir1(XoverOrder,fcNorm,'low',chebwin((XoverOrder+1),ripple1), 'noscale');
filtHi = fir1(XoverOrder,fcNorm,'high',chebwin((XoverOrder+1),ripple1), 'noscale');
end

