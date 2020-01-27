disp('Question 2 is running ...');
%% global variables or settings
load('eegdata.mat')
fs=250;

t = ((1:2500)-1)/fs;
% task = data{1}{4};
% plot(t,task(1,:))

c3=1;c4=2;p3=3;p4=4;o1=5;o2=6;eog=7;


task_group = cell(size(data));
subject_group = cell(size(data));
trial_group = cell(size(data));

for i=1:length(data)
    subject_group{i}=data{1,i}{1,1};
    task_group{i}=data{1,i}{1,2};
    trial_group{i}=data{1,i}{1,3};
end

multip_task_idx = [];
rotation_task_idx = [];
for i=1:length(data)
    if isequal(task_group{i},'multiplication')
        multip_task_idx = [multip_task_idx i];
    end
    if isequal(task_group{i},'rotation')
        rotation_task_idx = [rotation_task_idx i];
    end
end

for i=1:length(multip_task_idx)
    multiplication_task_data{i}=data{1,multip_task_idx(i)}{1,4};
end
for i=1:length(rotation_task_idx)
    rotation_task_data{i}=data{1,rotation_task_idx(i)}{1,4};
end

DataNumRows = size(multiplication_task_data{1},1);
% Removing Noise from Multiplication Task Data
for i=1:length(multip_task_idx)
    data=multiplication_task_data{i};
    for j=1:DataNumRows
        [num,den] = iirnotch(60/(fs/2),0.25/(fs/2));
        temp(j,:) = filter(num,den,data(j,:));
    end
    multip_task_filtered{i}= temp;
end
figure;
subplot(211); pwelch(multiplication_task_data{1}(c3,:))
title('Multiplication Task : before Noise Cancellation (iirnotch)')
subplot(212); pwelch(multip_task_filtered{1}(c3,:))
title('Multiplication Task : before Noise Cancellation (iirnotch)')

% Removing Noise from Rotation Task Data
for i=1:length(rotation_task_idx)
    data=rotation_task_data{i};
    for j=1:DataNumRows
        [num,den]=iirnotch(60/(fs/2),0.25/(fs/2));
        temp(j,:) = filter(num,den,data(j,:));
    end
    rotation_task_filtered{i}= temp;
end
figure
subplot(211); pwelch(rotation_task_data{1}(c3,:))
title('Rotation Task : before Noise Cancellation (iirnotch)')
subplot(212); pwelch(rotation_task_filtered{1}(c3,:))
title('Rotation Task : after Noise Cancellation (iirnotch)')

% Removing EOG Noise from Multiplication Task Data
rlsFilt = dsp.RLSFilter(25);

for i=1:length(multip_task_idx)
    data = multip_task_filtered{i};
    for j=1:DataNumRows-1
        [y,e] = rlsFilt(data(eog,:)',(data(eog,:))');
        temp(j,:) = data(j,:) - y';
    end
    multip_task_filtered2{i}= temp;
end
figure
subplot(211); plot(multip_task_filtered{1}(c3,:))
title('Multiplication Task : before EOG Noise Cancellation')
subplot(212); plot(multip_task_filtered2{1}(c3,:))
title('Multiplication Task : after EOG Noise Cancellation')

% Removing EOG Noise from Rotation Task Data
for i=1:length(rotation_task_idx)
    data=rotation_task_filtered{i};
    for j=1:DataNumRows-1
        [y,e] = rlsFilt(data(eog,:)',(data(eog,:))');
        temp(j,:) = data(j,:)-y';
    end
    rotation_task_filtered2{i}= temp;
end
figure
subplot(211); plot(rotation_task_filtered{1}(1,:))
title('Rotation Task Data before EOG Noise Cancellation')
subplot(212); plot(rotation_task_filtered2{1}(1,:))
title('Rotation Task Data after EOG Noise Cancellation')

%% Dividing Dataset to Train and Test subsets
% idx = randperm(z);
% TrainX = (X(idx(1:round(Ptrain.*z)),:))';
% TrainY = (Y(idx(1:round(Ptrain.*z)),:))';
% TestX = (X(idx(round(Ptrain.*z)+1:end),:))';
% TestY = (Y(idx(round(Ptrain.*z)+1:end),:))';

classes =[zeros(1,length(multip_task_idx)) ones(1,length(rotation_task_idx))];
rand_idx=randperm(2*length(multip_task_idx));
train_ind=rand_idx(1:floor(0.8*length(rand_idx)));
train_classess=classes(train_ind);

test_ind=rand_idx(floor(0.8*length(rand_idx))+1:end);
test_classess=classes(test_ind);
train_data_all=[multip_task_filtered2 rotation_task_filtered2];

for i=1:length(train_ind)
    train_data{i}=train_data_all{(train_ind(i))}; 
end
for i=1:length(test_ind)
    test_data{i}=train_data_all{(test_ind(i))}; 
end


%% A, B
delta_feature = zeros(104,6);
theta_feature = zeros(104,6);
alpha_feature = zeros(104,6);
beta_feature  = zeros(104,6);
f = 125/2*linspace(0,1,1024/2+1);
delta_f = f(2)-f(1);
ch = [3, 5, 8, 12, 19, 21];
for i = 1:2
    for j = 1:86
        for k = 1:6
        pxx = pwelch(class(ch(k),625:1000,86*(i-1)+j),125,100,1024);
        delta_feature(86*(i-1)+j,k) = sum(pxx(1:32).^2)*delta_f;
        theta_feature(86*(i-1)+j,k) = sum(pxx(33:66).^2)*delta_f;
        alpha_feature(86*(i-1)+j,k) = sum(pxx(67:107).^2)*delta_f;
        beta_feature(86*(i-1)+j,k)  = sum(pxx(108:247).^2)*delta_f;
        end
    end
end
%%%%%% Asymmetry_ratio
asymmetry_ratio_delta_feature = zeros(104,9);
asymmetry_ratio_theta_feature = zeros(104,9);
asymmetry_ratio_alpha_feature = zeros(104,9);
asymmetry_ratio_beta_feature  = zeros(104,9);
pair_ch = [3 5;3 12;3 21;8 5;8 12;8 21;19 5;19 12;19 21];
for i = 1:2
for j = 1:86
for k = 1:9
ch_R = find(ch == pair_ch(k,2));
ch_L = find(ch == pair_ch(k,1));
asymmetry_ratio_delta_feature(86*(i-1)+j,k) = (delta_feature(86*(i-1)+j,ch_R)-delta_feature(86*(i-1)+j,ch_L))/(delta_feature(86*(i-1)+j,ch_R)+delta_feature(86*(i-1)+j,ch_L));
asymmetry_ratio_theta_feature(86*(i-1)+j,k) = (theta_feature(86*(i-1)+j,ch_R)-theta_feature(86*(i-1)+j,ch_L))/(theta_feature(86*(i-1)+j,ch_R)+theta_feature(86*(i-1)+j,ch_L));
asymmetry_ratio_alpha_feature(86*(i-1)+j,k) = (alpha_feature(86*(i-1)+j,ch_R)-alpha_feature(86*(i-1)+j,ch_L))/(alpha_feature(86*(i-1)+j,ch_R)+alpha_feature(86*(i-1)+j,ch_L));
asymmetry_ratio_beta_feature(86*(i-1)+j,k)  = (beta_feature(86*(i-1)+j,ch_R)-beta_feature(86*(i-1)+j,ch_L))/(beta_feature(86*(i-1)+j,ch_R)+beta_feature(86*(i-1)+j,ch_L));
end
end
end
all_features = [delta_feature,theta_feature,alpha_feature,beta_feature,asymmetry_ratio_delta_feature,asymmetry_ratio_theta_feature,asymmetry_ratio_alpha_feature,asymmetry_ratio_beta_feature];
v=1:1:104;
v=v(randperm(104));
all_features_shuffle = all_features(v,:);
all_features_train   = all_features(1:115,:);
all_features_test    = all_features(116:104,:);
label= ones(104,1);
label(73:104)  = 2*ones(86,1);
label = label(v,:);
label_train = label(1:115,:);
label_test  = label(116:104);
%%  C
index_class_1 = find(label_train==1);
index_class_2 = find(label_train==2);
m1 = mean(all_features_train(index_class_1,:),1);
m2 = mean(all_features_train(index_class_2,:),1);
var1 = var(all_features_train(index_class_1,:),1);
var2 = var(all_features_train(index_class_2,:),1);
m = mean(all_features_train,1);
fisher_score = zeros(1,60);
for i = 1:60
   fisher_score(1,i) =  [(m1(i)-m(i))^2 + (m2(i)-m(i))^2]/(var1(i)+var2(i));
end
first_feature = find(fisher_score==max(fisher_score));
index = 1:60;
index(first_feature) = [];
fisher_score = zeros(1,59);
for j = 1:59
    a = [all_features_train(:,first_feature),all_features_train(:,index(j))];
    m1 = mean(a(index_class_1,:),1);
    m2 = mean(a(index_class_2,:),1);
    var1 = var(a(index_class_1,:),1);
    var2 = var(a(index_class_2,:),1);
    m = mean(a,1);
    S_b = (m1' - m')*(m1' - m')' + (m2' - m')*(m2' - m')';
    S1 = zeros(2,2);
    S2 = zeros(2,2);
    for i = 1:length(index_class_1);
        S1 = S1 + (a(index_class_1(i),:)' - m1')*(a(index_class_1(i),:)' - m1')';
    end
    for i = 1:length(index_class_2)
        S2 = S2 + (a(index_class_2(i),:)' - m2')*(a(index_class_2(i),:)' - m2')';
    end
    % scattering within class
    S_w = [length(index_class_1)*S1 + length(index_class_2)*S2] / [length(index_class_1)+length(index_class_2)];
    fisher_score(j) = trace(S_b)/trace(S_w);
end
second_feature = index(find(fisher_score==max(fisher_score)));

index = 1:60;
index(second_feature) = [];
index(first_feature)  = [];
fisher_score = zeros(1,58);
for j = 1:58
    a = [all_features_train(:,first_feature),all_features_train(:,second_feature),all_features_train(:,index(j))];
    m1 = mean(a(index_class_1,:),1);
    m2 = mean(a(index_class_2,:),1);
    var1 = var(a(index_class_1,:),1);
    var2 = var(a(index_class_2,:),1);
    m = mean(a,1);
    S_b = (m1' - m')*(m1' - m')' + (m2' - m')*(m2' - m')';
    S1 = zeros(3,3);
    S2 = zeros(3,3);
    for i = 1:length(index_class_1)
        S1 = S1 + (a(index_class_1(i),:)' - m1')*(a(index_class_1(i),:)' - m1')';
    end
    for i = 1:length(index_class_2)
        S2 = S2 + (a(index_class_2(i),:)' - m2')*(a(index_class_2(i),:)' - m2')';
    end
    % scattering within class
    S_w = [length(index_class_1)*S1 + length(index_class_2)*S2] / [length(index_class_1)+length(index_class_2)];
    fisher_score(j) = trace(S_b)/trace(S_w);
end
third_feature = index(find(fisher_score==max(fisher_score)));

%%  D
syms b1
syms b2
syms b3
cov_class_1 = cov(all_features_train(index_class_1,[first_feature, second_feature, third_feature]));
cov_class_2 = cov(all_features_train(index_class_2,[first_feature, second_feature, third_feature]));

m1 = mean(all_features_train(index_class_1,[first_feature, second_feature, third_feature]),1)';
m2 = mean(all_features_train(index_class_2,[first_feature, second_feature, third_feature]),1)';

% boundary of first class
A_1 = (-0.5)*inv(cov_class_1);
b_1 = inv(cov_class_1) * m1;
c_1 = (-0.5)* m1' * inv(cov_class_1) * m1 - 0.5 * log(det(cov_class_1));
d_1 = [b1 b2 b3]*A_1*[b1; b2; b3] + b_1'*[b1; b2; b3] + c_1;
% boundary of second class
A_2 = (-0.5)*inv(cov_class_2);
b_2 = inv(cov_class_2) * m2;
c_2 = (-0.5)* m2' * inv(cov_class_2) * m2 - 0.5 * log(det(cov_class_2));
d_2 = [b1 b2 b3]*A_2*[b1; b2; b3] + b_2'*[b1; b2; b3] + c_2;

figure();
plot3(all_features_train(index_class_1,first_feature),all_features_train(index_class_1,second_feature),all_features_train(index_class_1,third_feature),'r*');
hold on
plot3(all_features_train(index_class_2,first_feature),all_features_train(index_class_2,second_feature),all_features_train(index_class_2,third_feature),'g*');
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('class 1','class 2')

% bayes classification
pred_bayes = zeros(1,29);
for i = 1:29
    bayes_value = zeros(1,2);
    bayes_value(1) = subs(d_1, [b1, b2, b3], [all_features_test(i,first_feature), all_features_test(i,second_feature), all_features_test(i,third_feature)]);
    bayes_value(2) = subs(d_2, [b1, b2, b3], [all_features_test(i,first_feature), all_features_test(i,second_feature), all_features_test(i,third_feature)]);
    pred_bayes(1,i) = find(bayes_value==max(bayes_value));
end
acc_bayes_with_three_features = sum(eq(label_test',pred_bayes))/29;
%  Euclidean distance classification
pred_Euclidean = zeros(1,29);
for i = 1:29
   Euclidean_value = zeros(1,2);
   Euclidean_value(1) = sqrt(sum((all_features_test(i,[first_feature,second_feature,third_feature])-m1').^2));
   Euclidean_value(2) = sqrt(sum((all_features_test(i,[first_feature,second_feature,third_feature])-m2').^2));
   pred_Euclidean(1,i) = find(Euclidean_value==min(Euclidean_value));
end
acc_Euclidean_with_three_features = sum(eq(label_test',pred_Euclidean))/29;
%  Mahalanobis distance classification
d_M_1 = -0.5*([b1;b2;b3]-m1)' * inv(cov_class_1) * ([b1;b2;b3]-m1);
d_M_2 = -0.5*([b1;b2;b3]-m2)' * inv(cov_class_2) * ([b1;b2;b3]-m2);
pred_Mahalanobis = zeros(1,29);
for i = 1:29
    Mahalanobis_value = zeros(1,2);
    Mahalanobis_value(1) = subs(d_M_1, [b1, b2, b3], [all_features_test(i,first_feature), all_features_test(i,second_feature), all_features_test(i,third_feature)]);
    Mahalanobis_value(2) = subs(d_M_2, [b1, b2, b3], [all_features_test(i,first_feature), all_features_test(i,second_feature), all_features_test(i,third_feature)]);
    pred_Mahalanobis(1,i) = find(Mahalanobis_value==min(Mahalanobis_value));
end
acc_Mahalanobis_with_three_features = sum(eq(label_test',pred_Mahalanobis))/29;
%%  E
covariance = cov(all_features_train(:,[first_feature,second_feature,third_feature]));
[eigen_vector ~] = eig(covariance);
W = eigen_vector(:,[1,2]);
new_data = W' * all_features_train(:,[first_feature,second_feature,third_feature])';
new_data_test = W' * all_features_test(:,[first_feature,second_feature,third_feature])';
figure();
plot(new_data(1,index_class_1),new_data(2,index_class_1),'r*');
hold on
plot(new_data(1,index_class_2),new_data(2,index_class_2),'g*');

cov_class_1 = cov(new_data(:,index_class_1)');
cov_class_2 = cov(new_data(:,index_class_2)');

m1 = mean(new_data(:,index_class_1)',1)';
m2 = mean(new_data(:,index_class_2)',1)';

% boundary of first class
A_1 = (-0.5)*inv(cov_class_1);
b_1 = inv(cov_class_1) * m1;
c_1 = (-0.5)* m1' * inv(cov_class_1) * m1 - 0.5 * log(det(cov_class_1));
d_1 = A_1(1,1)*b1^2 + A_1(2,2)*b2^2 + (A_1(1,2)+A_1(2,1))*b1*b2 + b_1(1,1)*b1 + b_1(2,1)*b2 + c_1;
% boundary of second class
A_2 = (-0.5)*inv(cov_class_2);
b_2 = inv(cov_class_2) * m2;
c_2 = (-0.5)* m2' * inv(cov_class_2) * m2 - 0.5 * log(det(cov_class_2));
d_2 = A_2(1,1)*b1^2 + A_2(2,2)*b2^2 + (A_2(1,2)+A_2(2,1))*b1*b2 + b_2(1,1)*b1 + b_2(2,1)*b2 + c_2;

%  plot boundary
hold on
ezplot(d_1==d_2,[1.5*min(new_data(1,:)),1.5*max(new_data(1,:)),1.5*min(new_data(2,:)),1.5*max(new_data(2,:))]);

title('Entropy')
legend('class 1','class 2')
% bayes classification
pred_bayes = zeros(1,29);
for i = 1:29
    bayes_value = zeros(1,2);
    bayes_value(1) = subs(d_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    bayes_value(2) = subs(d_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    pred_bayes(1,i) = find(bayes_value==max(bayes_value));
end
acc_bayes_with_two_Entropy_features = sum(eq(label_test',pred_bayes))/29;
%  Euclidean distance classification
pred_Euclidean = zeros(1,29);
for i = 1:29
   Euclidean_value = zeros(1,2);
   Euclidean_value(1) = sqrt(sum((new_data_test(:,i)-m1).^2));
   Euclidean_value(2) = sqrt(sum((new_data_test(:,i)-m2).^2));
   pred_Euclidean(1,i) = find(Euclidean_value==min(Euclidean_value));
end
acc_Euclidean_with_two_Entropy_features = sum(eq(label_test',pred_Euclidean))/29;
%  Mahalanobis distance classification
d_M_1 = -0.5*([b1;b2]-m1)' * inv(cov_class_1) * ([b1;b2]-m1);
d_M_2 = -0.5*([b1;b2]-m2)' * inv(cov_class_2) * ([b1;b2]-m2);

pred_Mahalanobis = zeros(1,29);
for i = 1:29
    Mahalanobis_value = zeros(1,2);
    Mahalanobis_value(1) = subs(d_M_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    Mahalanobis_value(2) = subs(d_M_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    pred_Mahalanobis(1,i) = find(Mahalanobis_value==min(Mahalanobis_value));
end
acc_Mahalanobis_with_two_Entropy_features = sum(eq(label_test',pred_Mahalanobis))/29;
%%  F
m1 = mean(all_features_train(index_class_1,[first_feature, second_feature, third_feature]),1)';
m2 = mean(all_features_train(index_class_2,[first_feature, second_feature, third_feature]),1)';
m  = mean(all_features_train(:,[first_feature, second_feature, third_feature]),1)';
S1 = zeros(3,3);
for i = 1:length(index_class_1)
   S1 = S1 + (all_features_train(index_class_1(i),[first_feature, second_feature, third_feature])'-m1)*(all_features_train(index_class_1(i),[first_feature, second_feature, third_feature])'-m1)';
end
S2 = zeros(3,3);
for i = 1:length(index_class_2)
   S2 = S2 + (all_features_train(index_class_2(i),[first_feature, second_feature, third_feature])'-m2)*(all_features_train(index_class_2(i),[first_feature, second_feature, third_feature])'-m2)';
end

S_W = [length(index_class_1)*S1 + length(index_class_2)*S2] / [length(index_class_1)+length(index_class_2)];
S_B = (m1-m)*(m1-m)' + (m2-m)*(m2-m)';
[eigen_vectors, eigen_values] = eig(S_W,S_B);
W_GEVD = eigen_vectors(:,2:3);
new_data = W_GEVD' * all_features_train(:,[first_feature,second_feature,third_feature])';
new_data_test = W_GEVD' * all_features_test(:,[first_feature,second_feature,third_feature])';
figure();
plot(new_data(1,index_class_1),new_data(2,index_class_1),'r*');
hold on
plot(new_data(1,index_class_2),new_data(2,index_class_2),'g*');
%  plot boundary
hold on
ezplot(d_1==d_2,[1.5*min(new_data(1,:)),1.5*max(new_data(1,:)),1.5*min(new_data(2,:)),1.5*max(new_data(2,:))]);
title('Fisher')
legend('class 1','class 2')
cov_class_1 = cov(new_data(:,index_class_1)');
cov_class_2 = cov(new_data(:,index_class_2)');
m1 = mean(new_data(:,index_class_1)',1)';
m2 = mean(new_data(:,index_class_2)',1)';
% boundary of first class
A_1 = (-0.5)*inv(cov_class_1);
b_1 = inv(cov_class_1) * m1;
c_1 = (-0.5)* m1' * inv(cov_class_1) * m1 - 0.5 * log(det(cov_class_1));
d_1 = A_1(1,1)*b1^2 + A_1(2,2)*b2^2 + (A_1(1,2)+A_1(2,1))*b1*b2 + b_1(1,1)*b1 + b_1(2,1)*b2 + c_1;
% boundary of second class
A_2 = (-0.5)*inv(cov_class_2);
b_2 = inv(cov_class_2) * m2;
c_2 = (-0.5)* m2' * inv(cov_class_2) * m2 - 0.5 * log(det(cov_class_2));
d_2 = A_2(1,1)*b1^2 + A_2(2,2)*b2^2 + (A_2(1,2)+A_2(2,1))*b1*b2 + b_2(1,1)*b1 + b_2(2,1)*b2 + c_2;

pred_bayes = zeros(1,29);
for i = 1:29
    bayes_value = zeros(1,2);
    bayes_value(1) = subs(d_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    bayes_value(2) = subs(d_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);

    pred_bayes(1,i) = find(bayes_value==max(bayes_value));
end
acc_bayes_with_two_Fisher_features = sum(eq(label_test',pred_bayes))/29;
%  Euclidean distance classification
pred_Euclidean = zeros(1,29);
for i = 1:29
   Euclidean_value = zeros(1,2);
   Euclidean_value(1) = sqrt(sum((new_data_test(:,i)-m1).^2));
   Euclidean_value(2) = sqrt(sum((new_data_test(:,i)-m2).^2));

   s = find(Euclidean_value==min(Euclidean_value));
   pred_Euclidean(1,i) = s(1);
end
acc_Euclidean_with_two_Fisher_features = sum(eq(label_test',pred_Euclidean))/29;
%  Mahalanobis distance classification
d_M_1 = -0.5*([b1;b2]-m1)' * inv(cov_class_1) * ([b1;b2]-m1);
d_M_2 = -0.5*([b1;b2]-m2)' * inv(cov_class_2) * ([b1;b2]-m2);

pred_Mahalanobis = zeros(1,29);
for i = 1:29
    Mahalanobis_value = zeros(1,2);
    Mahalanobis_value(1) = subs(d_M_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    Mahalanobis_value(2) = subs(d_M_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);

    pred_Mahalanobis(1,i) = find(Mahalanobis_value==min(Mahalanobis_value));
end
acc_Mahalanobis_with_two_Fisher_features = sum(eq(label_test',pred_Mahalanobis))/29;
%% G
covariance = cov(all_features_train(:,[first_feature,second_feature,third_feature]));
[eigen_vector, eigen_value] = eig(covariance);
W = eigen_vector(:,2:3);
new_data = W' * all_features_train(:,[first_feature,second_feature,third_feature])';
new_data_test = W' * all_features_test(:,[first_feature,second_feature,third_feature])';
figure();
plot(new_data(1,index_class_1),new_data(2,index_class_1),'r*');
hold on
plot(new_data(1,index_class_2),new_data(2,index_class_2),'g*');
%  plot boundary
hold on
ezplot(d_1==d_2,[1.5*min(new_data(1,:)),1.5*max(new_data(1,:)),1.5*min(new_data(2,:)),1.5*max(new_data(2,:))]);
title('PCA')
legend('class 1','class 2')
cov_class_1 = cov(new_data(:,index_class_1)');
cov_class_2 = cov(new_data(:,index_class_2)');
m1 = mean(new_data(:,index_class_1)',1)';
m2 = mean(new_data(:,index_class_2)',1)';
% boundary of first class
A_1 = (-0.5)*inv(cov_class_1);
b_1 = inv(cov_class_1) * m1;
c_1 = (-0.5)* m1' * inv(cov_class_1) * m1 - 0.5 * log(det(cov_class_1));
d_1 = A_1(1,1)*b1^2 + A_1(2,2)*b2^2 + (A_1(1,2)+A_1(2,1))*b1*b2 + b_1(1,1)*b1 + b_1(2,1)*b2 + c_1;
% boundary of second class
A_2 = (-0.5)*inv(cov_class_2);
b_2 = inv(cov_class_2) * m2;
c_2 = (-0.5)* m2' * inv(cov_class_2) * m2 - 0.5 * log(det(cov_class_2));
d_2 = A_2(1,1)*b1^2 + A_2(2,2)*b2^2 + (A_2(1,2)+A_2(2,1))*b1*b2 + b_2(1,1)*b1 + b_2(2,1)*b2 + c_2;
pred_bayes = zeros(1,29);
for i = 1:29
    bayes_value = zeros(1,2);
    bayes_value(1) = subs(d_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    bayes_value(2) = subs(d_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);

    pred_bayes(1,i) = find(bayes_value==max(bayes_value));
end
acc_bayes_with_two_PCA_features = sum(eq(label_test',pred_bayes))/29;
%  Euclidean distance classification
pred_Euclidean = zeros(1,29);
for i = 1:29
   Euclidean_value = zeros(1,2);
   Euclidean_value(1) = sqrt(sum((new_data_test(:,i)-m1).^2));
   Euclidean_value(2) = sqrt(sum((new_data_test(:,i)-m2).^2));

   pred_Euclidean(1,i) = find(Euclidean_value==min(Euclidean_value));
end
acc_Euclidean_with_two_PCA_features = sum(eq(label_test',pred_Euclidean))/29;
%  Mahalanobis distance classification
d_M_1 = -0.5*([b1;b2]-m1)' * inv(cov_class_1) * ([b1;b2]-m1);
d_M_2 = -0.5*([b1;b2]-m2)' * inv(cov_class_2) * ([b1;b2]-m2);
pred_Mahalanobis = zeros(1,29);
for i = 1:29
    Mahalanobis_value = zeros(1,2);
    Mahalanobis_value(1) = subs(d_M_1, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);
    Mahalanobis_value(2) = subs(d_M_2, [b1, b2], [new_data_test(1,i), new_data_test(2,i)]);

    pred_Mahalanobis(1,i) = find(Mahalanobis_value==min(Mahalanobis_value));
end
acc_Mahalanobis_with_two_PCA_features = sum(eq(label_test',pred_Mahalanobis))/29;