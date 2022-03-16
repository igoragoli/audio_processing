% Quatrième Partie du TL - Méthodes de Pursuit

%% Initial definitions

[f, fs] = audioread('société.wav');

% Make signal mono, crop at 2.5s 
f = mean(f, 2);  
f = f(1:round(2.5*fs)); 

% Resample
fs_desired = 32000;
[p, q] = rat(fs_desired/fs);
f = resample(f, p, q);
fs = fs_desired;
