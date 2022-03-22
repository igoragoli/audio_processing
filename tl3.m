% Troisième partie du TL
% Apprentissage supervisioné

close all; clear; clc;

%% Loading data

train_file = matfile("train_tl3.mat");
train_arr = train_file.train_tl3;

o1_f1_train = train_arr(:, 1);
o1_f2_train = train_arr(:, 2);
o2_f1_train = train_arr(:, 3);
o2_f2_train = train_arr(:, 4);
o3_f1_train = train_arr(:, 5);
o3_f2_train = train_arr(:, 6);

test_file = matfile("test_tl3.mat");
test_arr = test_file.test_tl3;

o1_f1_test = test_arr(:, 1);
o1_f2_test = test_arr(:, 2);
o2_f1_test = test_arr(:, 3);
o2_f2_test = test_arr(:, 4);
o3_f1_test = test_arr(:, 5);
o3_f2_test = test_arr(:, 6);

val_file = matfile("val_tl3.mat");
val_arr = val_file.val_tl3;

o1_f1_val = val_arr(:, 1);
o1_f2_val = val_arr(:, 2);
o2_f1_val = val_arr(:, 3);
o2_f2_val = val_arr(:, 4);
o3_f1_val = val_arr(:, 5);
o3_f2_val = val_arr(:, 6);

%% Training - Class 1 versus Class 2 + Class 3

sig2_1=0.3;
nu_1=0.01;

n1 = length(o1_f1_train); n2 = length(o2_f1_train); n3 = length(o3_f1_train);

disp('-----> Class 1 versus Class 2 + Class 3');

data1 = [[o1_f1_train; o2_f1_train; o3_f1_train] [o1_f2_train; o2_f2_train; o3_f2_train]]'; 
classes1 = [ones(1,n1) -ones(1,n2+n3)];

[alphaloqo1, yloqo1, bbbopt1] = uncontretousoutil(data1, n1, n2+n3, classes1, sig2_1, nu_1, false);
tmpv2_1 = alphaloqo1.*classes1';

%% Training - Class 2 versus Class 1 + Class 3

sig2_2=0.3;
nu_2=0.01;

disp('-----> Class 2 versus Class 1 + Class 3');

data2 = [[o2_f1_train; o1_f1_train; o3_f1_train] [o2_f2_train; o1_f2_train; o3_f2_train]]'; 
classes2 = [ones(1,n2) -ones(1,n1+n3)];

[alphaloqo2, yloqo2, bbbopt2] = uncontretousoutil(data2, n2, n1+n3, classes2, sig2_2, nu_2, false);
tmpv2_2 = alphaloqo2.*classes2';

%% Training - Class 3 versus Class 1 + Class 2

sig2_3=0.3;
nu_3=0.01;

disp('-----> Class 3 versus Class 1 + Class 2');

data3 = [[o3_f1_train; o1_f1_train; o2_f1_train] [o3_f2_train; o1_f2_train; o2_f2_train]]'; 
classes3 = [ones(1,n3) -ones(1,n1+n2)];

[alphaloqo3, yloqo3, bbbopt3] = uncontretousoutil(data3, n3, n1+n2, classes3, sig2_3, nu_3, false);
tmpv2_3 = alphaloqo3.*classes3';


%% Training error (training set)

feat_train = {[o1_f1_train o1_f2_train], [o2_f1_train o2_f2_train], [o3_f1_train o3_f2_train]};

success_count=0;
total_count=0;
for i=1:length(feat_train)

    feat = feat_train{i};
    total_count=total_count+length(feat);

    for j=1:length(feat)

        sample=feat(j, :);
        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample', data1, sig2_1)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample', data2, sig2_2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample', data3, sig2_3)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

training_error=(total_count-success_count)/total_count*100;
disp(" ")
disp("Training Error: " + training_error + "%")


%% Generalization error (test set)

feat_test = {[o1_f1_test o1_f2_test], [o2_f1_test o2_f2_test], [o3_f1_test o3_f2_test]};

success_count=0;
total_count=0;
for i=1:length(feat_test)

    feat = feat_test{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j, :);
        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample', data1, sig2_1)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample', data2, sig2_2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample', data3, sig2_3)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

generalization_error=(total_count-success_count)/total_count*100;
disp(" ")
disp("Generalization Error: " + generalization_error + "%")

%% Cross-validation error (validation set)

feat_val = {[o1_f1_val o1_f2_val], [o2_f1_val o2_f2_val], [o3_f1_val o3_f2_val]};

success_count=0;
total_count=0;
for i=1:length(feat_val)

    feat = feat_val{i};
    total_count=total_count+length(feat);
    for j=1:length(feat)

        sample=feat(j, :);
        clnew1 = tmpv2_1'*lagis_rbf_gaussien(sample', data1, sig2_1)' + bbbopt1;
        clnew2 = tmpv2_2'*lagis_rbf_gaussien(sample', data2, sig2_2)' + bbbopt2;
        clnew3 = tmpv2_3'*lagis_rbf_gaussien(sample', data3, sig2_3)' + bbbopt3;
        
        clnew = [clnew1, clnew2, clnew3];
        [~, max_idx] = max(clnew);

        if max_idx==i
            success_count=success_count+1;
        end

    end

end

cross_validation_error=(total_count-success_count)/total_count*100;
disp(" ")
disp("Cross-validation Error: " + cross_validation_error + "%")

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