% Troisième partie du TL
% Apprentissage supervisioné

close all; clear; clc;

%% Extracting features of all types of birds

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
      n = 2^20;
      %n = 2^nextpow2(L);
    
      Y = fft(X, n);
      P = abs(Y/L);
      P = P(1:n/2+1);
    
      [peaks, locs] = findpeaks(P, 'MinPeakHeight', 0.01);
      locs = locs*Fs/n;
        
      [loc, loc_idx] = min(locs);
      peak = peaks(loc_idx);
    
      feat_array(j, :) = [loc peak];
    
    end

features{i, 1} = feat_array;

end

%% Normalization

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

%% Plotting features

figure(1);
clf;
grid on;
hold on;
plot(o1_f1,o1_f2,'or');
plot(o2_f1,o2_f2,'*b');
plot(o3_f1,o3_f2,'+g');
hold off;
