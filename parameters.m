% Parameters for all the scripts.
clear;

% 'fft': Execution parameters for get_frequencies_fft()
%   'T_frame': How much time (s) is represented by each frame.
%   'T_step': How much time (s) is skipped at each step. Overlapping
%       between frames is determined by choosing 'T_frame' and 'T_step'.
%   'N_fft': Transform length. Has to be big enough to allow for
%       zero-padding.
%   'error_threshold': Percentage of the reference frequency below which an
%       error will not be considered.
%   'peak_tolerable_distance': Minimum distance between peaks in frames.
%   'peak_detection_threshold': Percentage of the maximum frame peak below
%       which a peak will not be considered.
params.fft.T_frame = 40e-3;
params.fft.T_step = 5e-3;
params.fft.N_fft = 2^16;  % Big enough to allow for zero-padding
params.fft.error_threshold = 0.05;
params.fft.peak_tolerable_distance = round(params.fft.N_fft/200);
params.fft.peak_detection_threshold = 0.5;

save('params.mat');