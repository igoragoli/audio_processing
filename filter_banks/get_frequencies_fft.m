function f = get_frequencies_fft(x, h, fs, N_frame, step, N_fft, ...
    ignore_frames_before, ignore_frames_after, params)
    % Get fundamental frequencies for each frame of a given signal.
    % Args:
    %   x: signal
    %   h: window (it has to be of length N_frame)
    %   fs: sampling frequency of x
    %   N_frame: number of samples per frame
    %   step: number of samples skipped at each loop
    %   N_fft: FFT length (must be higher than N_frame for zero padding)
    %   ignore_frames_before: ignore frames before this value
    %   ignore_frames_after: ignore frames after this value
    %   params: parameter data in 'params.mat'
    
    N = length(x);  % Length of the signal
    number_of_steps = ceil(N/step);
    f = zeros(1, number_of_steps);  % Fundamental frequency vector
    
    reached_end = false;
    for m = 1:number_of_steps - 1
        if m > ignore_frames_before && m < ignore_frames_after && ~reached_end
            % Calculate FFT of signal frame, ignore negative frequencies
            index_ini = (m - 1)*step + 1; 
            index_end = index_ini + N_frame - 1;
            if index_end > N
                x_frame = x(index_ini:end);
                reached_end = true;
            else
                x_frame = x(index_ini:index_end);
            end
            x_frame = x_frame.*h(1:length(x_frame));
            X_frame = fft(x_frame, N_fft);
            X_frame = fftshift(X_frame);
            if mod(N_fft, 2) == 0
                X_frame = X_frame(N_fft/2 + 1:end);
                f_vals = fs*(0:N_fft/2 - 1)/N_fft;
            else
                X_frame = X_frame(ceil(N_fft/2):end);
                f_vals = fs*(0:N_fft/2 - 1)/N_fft;
            end  

            % Get relevant peaks in the FFT      
            [pks, locs] = findpeaks(abs(X_frame()), f_vals, 'MinPeakDistance', params.fft.peak_tolerable_distance);
            filt = pks > params.fft.peak_detection_threshold*max(pks);
            locs = locs(filt);
            f(m) = locs(1);

            %plot(f_vals, abs(X_frame()), 'b'); hold on;
            %scatter(locs(1), pks(1), 'r'); hold off;                 
        end
    end
end 