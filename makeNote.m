function [x] = makeNote(n, T, Fs)
%MAKENOTE creates a note signal with x(t) = cos(wt)exp(-3t/tau). where n is
%the n-th's note, tau = T is the note's length and Fs is the sampling frequency.

f0 = 440; % [Hz]

t   = linspace(0, T, T*Fs);
f   = 2^(n/12)*f0;
tau = T;
x   = cos(2*pi*f*t).*exp(-3*t/tau);
end

