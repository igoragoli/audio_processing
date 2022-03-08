% Deuxième Partie du TL - Ondelettes
% Débruitage

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
if true
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

[pk, loc] = findpeaks(denoising_ratio, thresholds, 'NPeaks', 1, ...
    'MinPeakProminence', 0.2);
best_ratio = pk;
best_threshold = loc;

fprintf('\nThresholding type: %c', thresholding_type);
fprintf('\nbest_ratio: %.4f \nbest_threshold: %.4f\n', best_ratio, best_threshold);

%% Plotting denoising ratio
figure();
plot(thresholds, denoising_ratio, 'b');
hold on; stem(loc, pk, 'r'); hold off;
title('Denoising ratio graph');
xlabel('Thresholds');
ylabel('Denoising ratio');
xlim([0, 1]);

%% Denoising
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


%% Comparison with universal thresholding

universal_threshold = sqrt(2*log(length(x)))*median(abs(d{1}))/0.6745;

c_denoise_uni= wthresh(c, thresholding_type, universal_threshold);
y_denoise_uni = waverec(c_denoise_uni, l, wname);
error_denoise_uni = y_denoise_uni - x;
energy_error_denoise_uni = sum(error_denoise_uni.^2);
SNR_denoise_uni = snr(x, error_denoise_uni);

fprintf('\nComparison with universal thresholding.\n');
fprintf('energy_error_denoise_uni: %.4f\n', energy_error_denoise_uni);
fprintf('SNR_denoise_uni: %.4f\n', SNR_denoise_uni);


%% Comparison with Matlab wavelet denoising method
y_denoise_mat = wden(y, 'rigrsure', 's', 'sln', level, wname);

error_denoise_mat = y_denoise_mat - x;
energy_error_denoise_mat = sum(error_denoise_mat.^2);
SNR_denoise_mat = snr(x, error_denoise_mat);

fprintf('\nComparison with rigrsure thresholding.\n');
fprintf('energy_error_denoise_mat: %.4f\n', energy_error_denoise_mat);
fprintf('SNR_denoise_mat: %.4f\n', SNR_denoise_mat);

%% Notes
% Using 'db1' wavelets:
% SNR with soft thresholding: 13.7 dB, SNR with hard thresholding:
% 13.5 dB
% Although the SNR is over 2 dB higher than in manually soft thresholding,
% qualitatively, the denoised sounds are very similar.

% Using 'sym20' wavelets:
% SNR with soft thresholding: 20.2881 dB, SNR with hard thresholding: 
% 23.9064 dB
% SNR with rigrsure thresholding: 24.2458 dB

% TODO : Suggestion from the prof, calculate things with other wavelets. 
% 'db1', 'sym20', 'sym18', and one other I think.