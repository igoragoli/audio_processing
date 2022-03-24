%% Initial definitions
close all; clear; clc;
fs = 16000;
Ts = 1/fs;
T = 0.2;
n = 0:1:round(T/Ts) - 1;

%% Original signal to be studied
f0 = 500; wn0 = 2*pi*f0/fs; phi0 = 0;
sig2_v = 0.01;

x = 0.5*sin(wn0*n + phi0);  % Pure signal
v = randn(size(x))*sqrt(sig2_v);  % Noise
y = x + v;  % Noisy signal

%% Plotting scalograms of x, v, t with Continuous Wavelet Transform
wname_cwt = 'morse';
if false
    figure();
    cwt(x, wname_cwt);
    title('Scalogram of x(t)');
    figure();
    cwt(v, wname_cwt);
    title('Scalogram of v(t)');
    figure();
    cwt(y, wname_cwt);
    title('Scalogram of y(t) = x(t) + v(t)');
end

%% Decomposing signal with Discrete Wavelet Transform
wname = 'sym20';
level = 10;
fprintf('Wavelet type: %s\n', wname);
fprintf('Number of filter bank levels: %d\n', level);
[c, l] = wavedec(y, level, wname);
d = detcoef(c, l, 1:1:level);

if false
    figure();
    level_plot = 5;
    subplot(level_plot + 1, 1, 1);
    plot(y);
    title('Original signal');
    for i = 1:level_plot
        subplot(level_plot + 1, 1, i + 1);
        plot(d{i}, 'b');
        title(sprintf('Level %d Detail Coefficients', i)); 
    end
end

%% Plotting coefficient histogram
coef_vec = [d{:}];
coef_vec = abs(coef_vec);
edges = 0:0.1:3;

figure();
histogram(coef_vec, edges, 'FaceColor', 'b');
xlim([0, 3]);
title('Detail coefficient histogram');

%% Wavelet thresholding analysis through denoising ratio

number_of_realizations = 1; % Number of realizations
thresholding_type = 'h'; 
thresholds = 0:0.01:1;
denoising_ratio = zeros(1, length(thresholds));
fprintf('\nPerforming wavelet thresholding analysis\n');
for l = 1:number_of_realizations
    fprintf('%d/%d\n', l, number_of_realizations);
    % Draw at each realization a new statistically independent random
    % number sequence
    v = randn(size(x))*sqrt(sig2_v);  % Noise
    y = x + v;  % Noisy signal 
    
    % Calculate discrete wavelet transform coefficients
    [c, l] = wavedec(y, level, wname);
    
    for i = 1:length(thresholds)
        % Threshold the coefficients through soft thresholding
        c_denoise = wthresh(c, thresholding_type, thresholds(i));

        % Reconstruct the denoised signal
        y_denoise = waverec(c_denoise, l, wname);

        % Calculate the denoising ratio for the current threshold
        error_noise = y - x;
        energy_error_noise = sum(error_noise.^2);
        error_denoise = y_denoise - x;
        energy_error_denoise = sum(error_denoise.^2);
        
        % Update denoising ratio mean
        denoising_ratio(i) = denoising_ratio(i) + ...
            (energy_error_noise/energy_error_denoise)/number_of_realizations;
    end
end

%% Plotting denoising ratio
figure();
plot(thresholds, denoising_ratio, 'b');
title('Denoising ratio graph');
xlabel('Thresholds');
ylabel('Denoising ratio');
xlim([0, 1]);

%% Denoising
best_threshold = 0.4;
c_denoise = wthresh(c, thresholding_type, best_threshold);
y_denoise = waverec(c_denoise, l, wname);

error_noise = y - x;
energy_error_noise = sum(error_noise.^2);
error_denoise = y_denoise - x;
energy_error_denoise = sum(error_denoise.^2);

SNR_noise = snr(x, error_noise);
SNR_denoise = snr(x, error_denoise);

fprintf('\nenergy_error_noise: %.4f \nenergy_error_denoise: %.4f\n', ...
    energy_error_noise, energy_error_denoise);
fprintf('SNR_noise: %.4f \nSNR_denoise: %.4f\n', SNR_noise, SNR_denoise);

%% Plotting denoised signal
figure();
subplot(2, 1, 1);
plot(n, x, 'g', n, y_denoise, 'b', n, y, 'r');
legend('Clean signal', 'Denoised signal', 'Noisy signal');
subplot(2, 1, 2);
plot(n, x, 'g', n, y_denoise, 'b', n, y, 'r');
legend('Clean signal', 'Denoised signal', 'Noisy signal');
xlim([0, 200]);

figure();
subplot(2, 1, 1);
plot(n, error_noise, 'r', n, error_denoise, 'b');
legend('Noisy sig. error', 'Denoised sig. error');
subplot(2, 1, 2);
plot(n, error_noise.^2, 'r', n, error_denoise.^2, 'b');
legend('Noisy sig. SE', 'Denoised sig. SE');

%% Plotting scalograms of the denoised signal with Continuous Wavelet Transform
wname_cwt = 'morse';
figure();
cwt(y_denoise, wname_cwt);
title('Scalogram of y_{denoise}(t)');
