% Troisième partie du TL
% Apprentissage supervisioné

close all; clear; clc;

%% Extracting features of all types of birds ("SONS")

dirs_oiseaux = dir("SONS");
dirs_oiseaux = dirs_oiseaux(~ismember({dirs_oiseaux.name},{'.','..'}));
num_dirs = length(dirs_oiseaux);

features = cell(num_dirs, 1);
for i = 1:length(dirs_oiseaux)

    dir_oiseau = dirs_oiseaux(i).name;
    files = dir(fullfile("SONS", dir_oiseau, '*.wav'));
    num_files = length(files);
    
    feat_array = zeros(num_files, 2);
    for j = 1:num_files
      filename = files(j).name;
      fullfilename = fullfile("SONS", dir_oiseau, filename);
    
      [X, Fs] = audioread(fullfilename);
      L = length(X);
      
      n = 2^20; % Number of points for the FFT
      Y = fft(X, n);
      P = abs(Y/L);
      P = P(1:n/2+1);
    
      % Finding the peaks of the FFT
      [peaks, locs] = findpeaks(P, 'MinPeakHeight', 0.01);
      locs = locs*Fs/n;
      [loc, loc_idx] = min(locs);
      peak = peaks(loc_idx);
    
      feat_array(j, :) = [loc peak];
    
    end

features{i, 1} = feat_array;

end

%% Normalization ("SONS")

feat1 = features{1, 1};
feat2 = features{2, 1};
feat3 = features{3, 1};

o1_f1 = feat1(:, 1); o1_f2 = feat1(:, 2);
o2_f1 = feat2(:, 1); o2_f2 = feat2(:, 2);
o3_f1 = feat3(:, 1); o3_f2 = feat3(:, 2);

mean1=mean([o1_f1; o2_f1; o3_f1]);
mean2=mean([o1_f2; o2_f2; o3_f2]);
std1=std([o1_f1; o2_f1; o3_f1]);
std2=std([o1_f2; o2_f2; o3_f2]);

o1_f1=(o1_f1-mean1)/std1;
o1_f2=(o1_f2-mean2)/std2;
o2_f1=(o2_f1-mean1)/std1;
o2_f2=(o2_f2-mean2)/std2;
o3_f1=(o3_f1-mean1)/std1;
o3_f2=(o3_f2-mean2)/std2;

%% Plotting features distribution ("SONS")

% figure(1);
% clf;
% grid on;
% hold on;
% plot(o1_f1,o1_f2,'or');
% plot(o2_f1,o2_f2,'*b');
% plot(o3_f1,o3_f2,'+g');
% title("Classes distribuition")
% xlabel("Normalized frequency")
% ylabel("Normalized amplitude")
% legend("Bird 1", "Bird 2", "Bird 3")
% hold off;

%% Division between training set and test set (80/20)

o1_f1_train = o1_f1(1:round(length(o1_f1)*0.8));
o1_f2_train = o1_f2(1:round(length(o1_f2)*0.8));
o2_f1_train = o2_f1(1:round(length(o2_f1)*0.8));
o2_f2_train = o2_f2(1:round(length(o2_f2)*0.8));
o3_f1_train = o3_f1(1:round(length(o3_f1)*0.8));
o3_f2_train = o3_f2(1:round(length(o3_f2)*0.8));

o1_f1_test = o1_f1(round(length(o1_f1)*0.2):end);
o1_f2_test = o1_f2(round(length(o1_f2)*0.2):end);
o2_f1_test = o2_f1(round(length(o2_f1)*0.2):end);
o2_f2_test = o2_f2(round(length(o2_f2)*0.2):end);
o3_f1_test = o3_f1(round(length(o3_f1)*0.2):end);
o3_f2_test = o3_f2(round(length(o3_f2)*0.2):end);

%% Saving arrays in .mat files

train_tl3 = [o1_f1_train, o1_f2_train, o2_f1_train, o2_f2_train, o3_f1_train, o3_f2_train];
test_tl3 = [o1_f1_test, o1_f2_test, o2_f1_test, o2_f2_test, o3_f1_test, o3_f2_test];

save("train_tl3.mat", "train_tl3");
save("test_tl3.mat", "test_tl3");

%% Extracting features of all types of birds ("SONS-VC")

dirs_oiseaux = dir("SONS-VC");
dirs_oiseaux = dirs_oiseaux(~ismember({dirs_oiseaux.name},{'.','..'}));
num_dirs = length(dirs_oiseaux);

features = cell(num_dirs, 1);
for i = 1:length(dirs_oiseaux)

    dir_oiseau = dirs_oiseaux(i).name;
    files = dir(fullfile("SONS-VC", dir_oiseau, '*.wav'));
    num_files = length(files);
    
    feat_array = zeros(num_files, 2);
    for j = 1:num_files
      filename = files(j).name;
      fullfilename = fullfile("SONS-VC", dir_oiseau, filename);
    
      [X, Fs] = audioread(fullfilename);
      L = length(X);

      n = 2^20; % Number of points for the FFT
      Y = fft(X, n);
      P = abs(Y/L);
      P = P(1:n/2+1);
    
      % Finding the peaks of the FFT
      [peaks, locs] = findpeaks(P, 'MinPeakHeight', 0.01);
      locs = locs*Fs/n;
      [loc, loc_idx] = min(locs);
      peak = peaks(loc_idx);
    
      feat_array(j, :) = [loc peak];
    
    end

features{i, 1} = feat_array;

end

%% Normalization ("SONS-VC")

feat1 = features{1, 1};
feat2 = features{2, 1};
feat3 = features{3, 1};

o1_f1_val = feat1(:, 1); o1_f2_val = feat1(:, 2);
o2_f1_val = feat2(:, 1); o2_f2_val = feat2(:, 2);
o3_f1_val = feat3(:, 1); o3_f2_val = feat3(:, 2);

mean1=mean([o1_f1_val; o2_f1_val; o3_f1_val]);
mean2=mean([o1_f2_val; o2_f2_val; o3_f2_val]);
std1=std([o1_f1_val; o2_f1_val; o3_f1_val]);
std2=std([o1_f2_val; o2_f2_val; o3_f2_val]);

o1_f1_val=(o1_f1_val-mean1)/std1;
o1_f2_val=(o1_f2_val-mean2)/std2;
o2_f1_val=(o2_f1_val-mean1)/std1;
o2_f2_val=(o2_f2_val-mean2)/std2;
o3_f1_val=(o3_f1_val-mean1)/std1;
o3_f2_val=(o3_f2_val-mean2)/std2;

%% Saving arrays in .mat files

val_tl3 = [o1_f1_val, o1_f2_val, o2_f1_val, o2_f2_val, o3_f1_val, o3_f2_val];

save("val_tl3.mat", "val_tl3");