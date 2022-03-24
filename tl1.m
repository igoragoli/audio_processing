% Première Partie du TL - Bancs de Filtres

%% Initial definitions
close all; clear; clc;
parameters;
include_vibrato = true;
Nfft = 2^12;
fs = 8000;
Ts = 1/fs;
T = 4;
n = 0:1:round(T/Ts) - 1;

%% Original signal to be studied
f0 = 880;
wn0 = 2*pi*f0/fs;
phi0 = [0 0 0 0];

% Configuring the vibrato
A_vib = 20;
f_vib = 5;
wn_vib = 2*pi*f_vib/fs;
phi_vib = 0; 
if include_vibrato
    vib = A_vib/f_vib*sin(wn_vib*n + phi_vib);
else
    vib = 0;
end

x = sin(1*wn0*n + phi0(1) + 1*vib) + ...
    sin(2*wn0*n + phi0(2) + 2*vib) + ...
    sin(3*wn0*n + phi0(3) + 3*vib) + ...
    sin(4*wn0*n + phi0(4) + 4*vib);

% Visualization of the Fourier Transform of the original signal
X = fft(x, Nfft);
X = fftshift(X);
f_vals = (-length(X)/2:length(X)/2 - 1)/length(X);
f_ticks = (-4:4) * f0;
figure();
plot(f_vals*fs, abs(X), 'k');
title('Original signal');
subtitle('|X(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

%% First decomposition
fc = fs/2;
wnc = 2*pi*fc/(2*pi*fs);  % Normalized cutoff frequency
xl = lowpass(x, wnc);
xh = highpass(x, wnc);

% Downsampling
N = 2;
fs_y = fs/N;
yl = xl(1:N:end); 
yh = xh(1:N:end);

% Visualization
Yl = fft(yl, Nfft);
Yl = fftshift(Yl);
f_vals = (-length(Yl)/2:length(Yl)/2 - 1)/length(Yl);
f_ticks = (-2:2) * f0;
figure();
subplot(2, 1, 1);
plot(f_vals*fs_y, abs(Yl), 'b');
title('First decomposition');
subtitle('|Y_l(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

Yh = fft(yh, Nfft);
Yh = fftshift(Yh);
f_vals = (-length(Yh)/2:length(Yh)/2 - 1)/length(Yh);
f_ticks = [480 1360]; f_ticks = [-flip(f_ticks), f_ticks];
subplot(2, 1, 2);
plot(f_vals*fs_y, abs(Yh), 'b');
subtitle('|Y_h(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

%% Second decomposition
fc = fs_y/N;
wnc = 2*pi*fc/(2*pi*fs_y);  % Normalized cutoff frequency
yll = lowpass(yl, wnc);
ylh = highpass(yl, wnc);
yhl = lowpass(yh, wnc);
yhh = highpass(yh, wnc);

% Downsampling
fs_y = fs_y/N;
yll = yll(1:N:end);
ylh = ylh(1:N:end);
yhl = yhl(1:N:end);
yhh = yhh(1:N:end);

% Visualization
Yll = fft(yll, Nfft);
Yll = fftshift(Yll);
f_vals = (-length(Yll)/2:length(Yll)/2 - 1)/length(Yll);
f_ticks = (-1:1) * f0;
figure();
subplot(4, 1, 1);
plot(f_vals*fs_y, abs(Yll), 'r');
title('Second decomposition');
subtitle('|Y_{ll}(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

Ylh = fft(ylh, Nfft);
Ylh = fftshift(Ylh);
f_vals = (-length(Ylh)/2:length(Ylh)/2 - 1)/length(Ylh);
f_ticks = [-240, 240];
subplot(4, 1, 2);
plot(f_vals*fs_y, abs(Ylh), 'r');
subtitle('|Y_{lh}(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

Yhl = fft(yhl, Nfft);
Yhl = fftshift(Yhl);
f_vals = (-length(Yhl)/2:length(Yhl)/2 - 1)/length(Yhl);
f_ticks = [-480, 480];
subplot(4, 1, 3);
plot(f_vals*fs_y, abs(Yhl), 'r');
subtitle('|Y_{hl}(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;

Yhh = fft(yhh, Nfft);
Yhh = fftshift(Yhh);
f_vals = (-length(Yhh)/2:length(Yhh)/2 - 1)/length(Yhl);
f_ticks = [-640, 640];
subplot(4, 1, 4);
plot(f_vals*fs_y, abs(Yhh), 'r');
subtitle('|Y_{hh}(e^{j \omega})|');
ax = gca;
ax.XTick = f_ticks;
fig = gcf;

%% Short-Term Fourier Transform (Transformée Fourier de Court Terme)

% Set up FFT parameters
load('params.mat');
Ts_y = 1/fs_y;
N_fft = params.fft.N_fft;
N_frame = round(params.fft.T_frame/Ts_y);
step = round(params.fft.T_step/Ts_y);
N = length(x);
number_of_steps = ceil(N/step);
ignore_frames_before = floor(0/T*number_of_steps);
ignore_frames_after = ceil(T/T*number_of_steps);

% Set up window
h = hann(N_frame)';

% Get frequencies through STFT
fll = get_frequencies_fft(yll, h, fs_y, N_frame, step, N_fft, ...
    ignore_frames_before, ignore_frames_after, params);
flh = get_frequencies_fft(ylh, h, fs_y, N_frame, step, N_fft, ...
    ignore_frames_before, ignore_frames_after, params);
fhl = get_frequencies_fft(yhl, h, fs_y, N_frame, step, N_fft, ...
    ignore_frames_before, ignore_frames_after, params);
fhh = get_frequencies_fft(yhh, h, fs_y, N_frame, step, N_fft, ...
    ignore_frames_before, ignore_frames_after, params);

% Plot graphs
T_step = T / length(fll);
t = 0:T_step:T - T_step;
figure();
set(0, 'DefaultAxesColorOrder', winter(4));
plot(t, fll, t, fhh, t, fhl, t, flh, 'LineWidth', 1);
legend('f_{ll}', 'f_{hh}', 'f_{hl}', 'f_{lh}');
title('f(t)');
