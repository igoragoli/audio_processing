% Troisième partie du TL
% Apprentissage supervisioné

close all; clear; clc;

%% Extracting features of all types of birds ("SONS")

dirs_oiseaux = dir("SONS");
dirs_oiseaux = dirs_oiseaux(~ismember({dirs_oiseaux.name},{'.','..'}));
%num_dirs = length(dirs_oiseaux);
num_dirs = round(length(dirs_oiseaux)/100);

features = cell(num_dirs, 1);
for i = 1:length(dirs_oiseaux)

    dir_oiseau = dirs_oiseaux(i).name;
    files = dir(fullfile("SONS", dir_oiseau, '*.wav'));
    %num_files = length(files);
    num_files = length(files)/20;
    
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

%% Plotting features

figure(1);
clf;
grid on;
hold on;
plot(o1_f1,o1_f2,'or');
plot(o2_f1,o2_f2,'*b');
plot(o3_f1,o3_f2,'+g');
hold off;

%% Division between training set and test set

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

%% Extracting features of all types of birds ("SONS-VC")

dirs_oiseaux = dir("SONS-VC");
dirs_oiseaux = dirs_oiseaux(~ismember({dirs_oiseaux.name},{'.','..'}));
num_dirs = length(dirs_oiseaux);

features = cell(num_dirs, 1);
for i = 1:length(dirs_oiseaux)

    dir_oiseau = dirs_oiseaux(i).name;
    files = dir(fullfile("SONS-VC", dir_oiseau, '*.wav'));
    %num_files = length(files);
    num_files = length(files)/20;
    
    feat_array = zeros(num_files, 2);
    for j = 1:num_files
      filename = files(j).name;
      fullfilename = fullfile("SONS-VC", dir_oiseau, filename);
    
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

%% Taining - Class 1 versus Class 2 + Class 3

sig2=0.5;
nu=0.01;
n1 = length(o1_f1_train); n2 = length(o2_f1_train); n3 = length(o3_f1_train);

disp('-----> Class 1 versus Class 2 + Class 3');

data1 = [[o1_f1_train; o2_f1_train; o3_f1_train] [o1_f2_train; o2_f2_train; o3_f2_train]]'; 
classes1 = [ones(1,n1) -ones(1,n2+n3)];

[alphaloqo1, yloqo1, bbbopt1] = uncontretousoutil(data1, n1, n2+n3, classes1, sig2, nu);
tmpv2_1 = alphaloqo1.*classes1';

%% Training - Class 2 versus Class 1 + Class 3

disp('-----> Class 2 versus Class 1 + Class 3');

data2 = [[o2_f1_train; o1_f1_train; o3_f1_train] [o2_f2_train; o1_f2_train; o3_f2_train]]'; 
classes2 = [ones(1,n2) -ones(1,n1+n3)];

[alphaloqo2, yloqo2, bbbopt2] = uncontretousoutil(data2, n2, n1+n3, classes2, sig2, nu);
tmpv2_2 = alphaloqo2.*classes2';

%% Training - Class 3 versus Class 1 + Class 2

disp('-----> Class 3 versus Class 1 + Class 2');

data3 = [[o3_f1_train; o1_f1_train; o2_f1_train] [o3_f2_train; o1_f2_train; o2_f2_train]]'; 
classes3 = [ones(1,n3) -ones(1,n1+n2)];

[alphaloqo3, yloqo3, bbbopt3] = uncontretousoutil(data3, n3, n1+n2, classes3, sig2, nu);
tmpv2_3 = alphaloqo3.*classes3';

%% Training error

feat_train = {[o1_f1_train o1_f2_train], [o2_f1_train o2_f2_train], [o3_f1_train o3_f2_train]};
feat_test = {[o1_f1_test o1_f2_test], [o2_f1_test o2_f2_test], [o3_f1_test o3_f2_test]};
feat_val = {[o1_f1_val o1_f2_val], [o2_f1_val o2_f2_val], [o3_f1_val o3_f2_val]};

success_count=0;
total_count=0;
for i=1:length(feat_train)

    feat = feat_train{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j);

        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample, data1, sig2)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample, data2, sig2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample, data3, sig2)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);
        
        disp("m " + max_idx + " i " + i)
        if max_idx==i
            success_count=success_count+1;
        end

    end

end

training_error=(total_count-success_count)/total_count;

%% Training error

feat_train = {[o1_f1_train o1_f2_train], [o2_f1_train o2_f2_train], [o3_f1_train o3_f2_train]};
feat_test = {[o1_f1_test o1_f2_test], [o2_f1_test o2_f2_test], [o3_f1_test o3_f2_test]};
feat_val = {[o1_f1_val o1_f2_val], [o2_f1_val o2_f2_val], [o3_f1_val o3_f2_val]};

success_count=0;
total_count=0;
for i=1:length(feat_train)

    feat = feat_train{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j);

        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample, data1, sig2)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample, data2, sig2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample, data3, sig2)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

training_error=(total_count-success_count)/total_count;

%% Generalization error

success_count=0;
total_count=0;
for i=1:length(feat_test)

    feat = feat_test{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j);

        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample, data1, sig2)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample, data2, sig2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample, data3, sig2)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

generalization_error=(total_count-success_count)/total_count;

%% Cross-validation error

success_count=0;
total_count=0;
for i=1:length(feat_val)

    feat = feat_val{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j);

        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample, data1, sig2)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample, data2, sig2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample, data3, sig2)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

cross_validation_error=(total_count-success_count)/total_count;

%% plot pour la route
% figure(7);
% clf;
% title('frontiere entre les 3 classes');
% hold on;
% imagesc(ppd1res,ppd2res,-clnew');
% grid on;
% hold on;
% plot(o1_f1,o1_f2,'*r');
% plot(o2_f1,o2_f2,'ob');
% plot(o3_f1,o3_f2,'+g');
% xlabel('dimension 1');
% ylabel('dimension 2');
% xlim([min(ppd1res) max(ppd1res)]);
% ylim([min(ppd2res) max(ppd2res)]);
% hold off;
% colormap('cool');