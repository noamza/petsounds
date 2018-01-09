function [s] = makeTune(n, T, Fs)
%MAKETUNE plays a tune composed by the vector of notes n and each note's
%length T
s = [];
for ii = 1:length(n)
    s = [s; (makeNote(n(ii), T(ii), Fs)).']; %#ok
end

sound(s, Fs);