% Deuxième Partie du TL - Ondelettes
% Séparation

%% Initial definitions
close all; clear; clc;
run('./ltfat-master/ltfatstart.m');  % Starts LTFAT toolkit

fs = 16000;
Ts = 1/fs;
T = 0.2;
n = 0:1:round(T/Ts) - 1;

%% Original signals to be studied
f0 = 100; wn0 = 2*pi*f0/fs; phi0 = 0; a1 = 0.5;
f1 = 500; wn1 = 2*pi*f1/fs; phi1 = 0; a2 = 0.5;
sig2_v = 0.01;

x0 = a1*sin(wn0*n + phi0);
x1 = a2*sin(wn1*n + phi1);
x = x0 + x1;

%% Decomposes signal with (Fast) Discrete Wavelet Transform
wname = 'sym20';
level = 10;
fprintf('\nWavelet type: %s\n', wname);
fprintf('Number of filter bank levels: %d\n', level);

% If true, uses Matlab's 'wavedec', otherwise, uses LTFAT's 'fwt'.
% Either way, informations from 'info' struct is taken from 'fwt' in order
% to use 'plotwavelets' in the following step.
use_wavedec = true;  

if use_wavedec
    % x
    [~, info] = fwt(x, wname, level);
    [c, l] = wavedec(x, level, wname);
    d = detcoef(c, l, 1:1:level);
    info.Lc = l(1:end-1)';
    
    % x0
    [~, info0] = fwt(x0, wname, level);
    [c0, l0] = wavedec(x0, level, wname);
    d0 = detcoef(c0, l, 1:1:level);
    info0.Lc = l0(1:end-1)';
    
    % x1
    [~, info1] = fwt(x1, wname, level);
    [c1, l1] = wavedec(x1, level, wname);
    d1 = detcoef(c, l1, 1:1:level);
    info1.Lc = l1(1:end-1)';
    
else
    % x
    [c, info] = fwt(x, wname, level);
    l = [info.Lc', info.Ls];
    d = detcoef(c, l, 1:1:level);
    
    % x0
    [c0, info0] = fwt(x0, wname, level);
    l0 = [info0.Lc', info0.Ls];
    d0 = detcoef(c0, l0, 1:1:level);
    
    % x1
    [c1, info1] = fwt(x1, wname, level);
    l1 = [info1.Lc', info1.Ls];
    d1 = detcoef(c1, l1, 1:1:level);
end

%% Plots scalogram of x[n]
figure();
plotwavelets(c, info, fs, 'dynrange', 100);
title('Scalogram of x[n]');
ax = gca; 
ax.FontSize = 9;

%% Plots detail coefficients for each level
if false
    figure();
    level_plot = 10;
    subplot(level_plot + 1, 1, 1);
    plot(x);
    title('Original signal');
    for i = 1:level_plot
        subplot(level_plot + 1, 1, i + 1);
        plot(d{i}, 'b');
        title(sprintf('Level %d Detail Coefficients', i)); 
    end
end

%% Splits detail coefficients into two groups
% Calculates power on each detail coefficient level
powers_d = zeros(length(d), 1);
for i = 1:length(d)
    powers_d(i) = sum(d{i}.^2)/length(d{i});
end

[~, idx] = maxk(powers_d, 2);
mean_idx = round(mean(idx));
fprintf("\nDetail coefficient levels with greatest power: %d and %d\n", idx);
fprintf("Splitting detail coefficient levels at index %d (higher levels include %d)\n", ...
    mean_idx, mean_idx);
idx1 = l((end - mean_idx)+1:end - 1);
idx2 = l(1:mean_idx);
split_idx = sum(idx2);
end_idx = sum(idx1) + sum(idx2);

% Sets to 0 detail coefficients which should be discarded
c0_sep = c; c1_sep = c;
c1_sep(1:split_idx) = 0;
c0_sep(split_idx + 1:end_idx) = 0;


%% Plots scalogram of x0[n], x1[n] and their separated signals
figure();
subplot(2, 2, 1);
plotwavelets(c0, info0, fs, 'dynrange', 100); ax = gca; ax.FontSize = 9;
title('Scalogram of x0[n]');

subplot(2, 2, 2);
plotwavelets(c1, info1, fs, 'dynrange', 100); ax = gca; ax.FontSize = 9;
title('Scalogram of x1[n]');

subplot(2, 2, 3);
plotwavelets(c0_sep, info, fs, 'dynrange', 100); ax = gca; ax.FontSize = 9;
title('Scalogram of x0_{sep}[n]');

subplot(2, 2, 4);
plotwavelets(c1_sep, info, fs, 'dynrange', 100); ax = gca; ax.FontSize = 9;
title('Scalogram of x1_{sep}[n]');


%% Reconstructs separated coefficients
x0_sep = waverec(c0_sep, l, wname);
x1_sep = waverec(c1_sep, l, wname);

error0 = x0_sep - x0;
error1 = x1_sep - x1;
SNR0= snr(x0, error0);
SNR1 = snr(x1, error1);

fprintf('\nSNR0: %.4f, SNR1: %.4f\n', SNR0, SNR1);

%% Plotting separated signals
figure();
subplot(2, 1, 1);
plot(n, x0, 'g', n, x0_sep, 'b');
title('Original and separated signals');
legend('x0[n]', 'x0_{sep}[n]');
xlim([200, 300]);
subplot(2, 1, 2);
plot(n, x1, 'g', n, x1_sep, 'b');
legend('x1[n]', 'x1_{sep}[n]');
xlim([200, 300]);

% Border effects
figure();
subplot(2, 1, 1);
plot(n, x0, 'g', n, x0_sep, 'b');
title('Original and separated signals (border effects)');
legend('x0[n]', 'x0_{sep}[n]');
xlim([0, 100]);
subplot(2, 1, 2);
plot(n, x1, 'g', n, x1_sep, 'b');
legend('x1[n]', 'x1_{sep}[n]');
xlim([0, 100]);

%% Plots errors
figure();
subplot(2, 1, 1);
title('Absolute errors');
plot(n, abs(error0), 'r');
xlim([-50, 3250]);
legend('|e_0[n]|');
subplot(2, 1, 2);
plot(n, abs(error1), 'r');
xlim([-50, 3250]);
legend('|e_1[n]|');
