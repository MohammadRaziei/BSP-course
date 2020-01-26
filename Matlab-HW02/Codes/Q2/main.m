disp('Question 2 is running ...');
%% global variables or settings
load('eegdata.mat')
fs=250;

t = ((1:2500)-1)/fs;
% task = data{1}{4};
% plot(t,task(1,:))

c3=1;c4=2;p3=3;p4=4;o1=5;o2=6;eog=7;


task_cell = cell(size(data));
subject_cell = cell(size(data));
trial_cell = cell(size(data));

for i=1:length(data)
    subject_cell{i}=data{1,i}{1,1};
    task_cell{i}=data{1,i}{1,2};
    trial_cell{i}=data{1,i}{1,3};
end

multip_task_idx = [];
rotation_task_idx = [];
for i=1:length(data)
    if isequal(task_cell{i},'multiplication')
        multip_task_idx = [multip_task_idx i];
    end
    if isequal(task_cell{i},'rotation')
        rotation_task_idx = [rotation_task_idx i];
    end
end

for i=1:length(multip_task_idx)
    multip_task_data{i}=data{1,multip_task_idx(i)}{1,4};
end
for i=1:length(rotation_task_idx)
    rotation_task_data{i}=data{1,rotation_task_idx(i)}{1,4};
end

DataNumRows = size(multip_task_data{1},1);
% Removing Noise from Multiplication Task Data
for i=1:length(multip_task_idx)
    data=multip_task_data{i};
    for j=1:DataNumRows
        [num,den]=iirnotch(60/(fs/2),0.25/(fs/2));
        temp(j,:) = filter(num,den,data(j,:));
    end
    multip_task_filtered{i}= temp;
end
figure;
subplot(211); pwelch(multip_task_data{1}(c3,:))
title('Multiplication Task Data before Noise Cancellation (iirnotch)')
subplot(212); pwelch(multip_task_filtered{1}(c3,:))
title('Multiplication Task Data before Noise Cancellation (iirnotch)')

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
title('Rotation Task Data before Noise Cancellation (iirnotch)')
subplot(212); pwelch(rotation_task_filtered{1}(c3,:))
title('Rotation Task Data after Noise Cancellation (iirnotch)')

% Removing EOG Noise from Multiplication Task Data
L=20;
rlsFilt = dsp.RLSFilter(L);

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
title('Multiplication Task Data before EOG Noise Cancellation')
subplot(212); plot(multip_task_filtered2{1}(c3,:))
title('Multiplication Task Data after EOG Noise Cancellation')

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

%% A
% first feature: Band Power
band_f0=[0 4 8 14];
band_f1=[3 7 13 20];
for i=1:floor(0.8*length(train_data_all))
    EEG=train_data{i};
    for j=1:6
        for k=1:4
            band_Power(i,j,k)=bandpower(EEG(j,:),fs,[band_f0(k) band_f1(k)]);
        end
    end
end
% second feature: Asymmetry Ratio
for i=1:0.8*length(train_data_all)
    c=1;
    for j=1:2:5
        for k=2:2:6
            for l=1:4
                R = band_Power(i,k,l);
                L = band_Power(i,j,l);
                AR(i,c)=(R-L)/(R+L);
                c=c+1;
            end
        end
    end
end
% Converting BW Power from rank 3 Tensor to rank 2 Tensor
for i=1:0.8*length(train_data_all)
    c=1;
    for j=1:6
        for k=1:4
            BW_Power_Feature(i,c)=band_Power(i,j,k);
            c=c+1;
        end
    end
end

% Features Combinations
all_features=[BW_Power_Feature AR];

%% B
% Feature selection
for i=1:60
    selected_feature=all_features(:,i);
    p_w0=length(find(train_classess==0));
    p_w1=length(find(train_classess==1));
    
    mu_0=mean(selected_feature(train_classess==0,:))';
    mu_1=mean(selected_feature(train_classess==1,:))';
    
    mu_total=p_w0*mu_0+p_w1*mu_1;
    %mu_total=mu_0+mu_1;
    mu_total=mu_total';
    Sb=p_w0*(mu_0-mu_total)*(mu_0-mu_total)' + (mu_1-mu_total)*(mu_1-mu_total)';
    x0=selected_feature(train_classess==0,:);
    x1=selected_feature(train_classess==1,:);
    cov0=0;
    for j=1:length(x0)
        %cov0=cov0+(x0(j,:)'-mu_0)*(x0(j,:)'-mu_0)';
        cov0=cov0+(x0(j,:)'-mu_total)*(x0(j,:)'-mu_total)';
    end
    cov1=0;
    for j=1:length(x1)
        %cov1=cov1+(x1(j,:)'-mu_1)*(x1(j,:)'-mu_1)';
        cov1=cov1+(x1(j,:)'-mu_total)*(x1(j,:)'-mu_total)';
    end

    %Sw=p_w0*cov0+p_w1*cov1;
    Sw=cov0+cov1;
    J(i)=trace(Sb)/trace(Sw);
end
best_feature=find(J==max(J));

% Finding the second best feature
k=1:60;
k=setdiff(k,best_feature);
J=0;
for i=k
    selected_feature=[all_features(:,best_feature) all_features(:,i) ];
    p_w0=length(find(train_classess==0));
    p_w1=length(find(train_classess==1));
    
    mu_0=mean(selected_feature(train_classess==0,:))';
    mu_1=mean(selected_feature(train_classess==1,:))';
    
    mu_total=p_w0*mu_0+p_w1*mu_1;
    %mu_total=mu_0+mu_1;
    mu_total=mu_total';
    Sb=p_w0*(mu_0-mu_total)*(mu_0-mu_total)' + (mu_1-mu_total)*(mu_1-mu_total)';
    x0=selected_feature(train_classess==0,:);
    x1=selected_feature(train_classess==1,:);
    % sum
    cov0=0;
    for j=1:length(x0)
        cov0=cov0+(x0(j,:)'-mu_total)*(x0(j,:)'-mu_total)';
    end
    % sum
    cov1=0;
    for j=1:length(x1)
        %cov1=cov1+(x1(j,:)'-mu_1)*(x1(j,:)'-mu_1)';
        cov1=cov1+(x1(j,:)'-mu_total)*(x1(j,:)'-mu_total)';
    end

    %Sw=p_w0*cov0+p_w1*cov1;
    Sw=cov0+cov1;
    J(i)=trace(Sb)/trace(Sw);
end
best_feature2=find(J==max(J));


% Finding the third best feature
k=1:60;
k=setdiff(k,[best_feature best_feature2]);
J=0;
for i=k
    selected_feature=[all_features(:,best_feature) all_features(:,best_feature2) all_features(:,i) ];
    p_w0=length(find(train_classess==0));
    p_w1=length(find(train_classess==1));
    
    mu_0=mean(selected_feature(train_classess==0,:))';
    mu_1=mean(selected_feature(train_classess==1,:))';
    
    mu_total=(p_w0*mu_0+p_w1*mu_1)';
    mu_total=mu_total';
    Sb=p_w0*(mu_0-mu_total)*(mu_0-mu_total)' + (mu_1-mu_total)*(mu_1-mu_total)';
    x0=selected_feature(train_classess==0,:);
    x1=selected_feature(train_classess==1,:);
    cov0=0;
    for j=1:length(x0)
        cov0=cov0+(x0(j,:)'-mu_total)*(x0(j,:)'-mu_total)';
    end
    cov1=0;
    for j=1:length(x1)
        cov1=cov1+(x1(j,:)'-mu_total)*(x1(j,:)'-mu_total)';
    end

    %Sw=p_w0*cov0+p_w1*cov1;
    Sw=cov0+cov1;
    J(i)=trace(Sb)/trace(Sw);
end
best_feature3=find(J==max(J));
all_best_features=[best_feature best_feature2 best_feature3];
reduced_data=[all_features(:,best_feature) all_features(:,best_feature2) all_features(:,best_feature3)];
figure
scatter3(reduced_data(train_classess==0,1),reduced_data(train_classess==0,2),reduced_data(train_classess==0,3))
hold on
scatter3(reduced_data(train_classess==1,1),reduced_data(train_classess==1,2),reduced_data(train_classess==1,3))

%% C
syms beta1 beta2 beta3
beta=[beta1 beta2 beta3]';
% Class 0
x0=reduced_data(train_classess==0,:);
mu0=mean(x0)';
sigma0=cov(x0);
A0=-0.5*inv(sigma0);
b0=inv(sigma0)*mu0;
c0=-0.5*mu0'*sigma0*mu0-0.5*det(sigma0);
d0=beta'*A0*beta + b0'*beta + c0;

% Class 1
x1=reduced_data(train_classess==1,:);
mu1=mean(x1)';
sigma1=cov(x1);
A1=-0.5*inv(sigma1);
b1=inv(sigma1)*mu1;
c1=-0.5*mu1'*sigma1*mu1-0.5*det(sigma1);
d1=beta'*A1*beta + b1'*beta + c1;
boundary=d0-d1;
boundary_bayes_bayesian=boundary;
boundary_euclidean_bayesian=sqrt(sum(beta-mu0).^2)-sqrt(sum(beta-mu1).^2);
boundary_mah_bayesian=-0.5*(beta-mu0)' * inv(sigma0) * (beta-mu0) + 0.5*(beta-mu1)' * inv(sigma1) * (beta-mu1);
figure
scatter3(reduced_data(train_classess==0,1),reduced_data(train_classess==0,2),reduced_data(train_classess==0,3))
hold on
scatter3(reduced_data(train_classess==1,1),reduced_data(train_classess==1,2),reduced_data(train_classess==1,3))
hold on
fimplicit3(boundary==0)
legend('Class 1','Class 2','Boundary')
%% D: Entropy
sigma=cov(reduced_data);
[V,lambda]=eig(sigma);
P=[V(:,1) V(:,2)];
P_entropy=P;
for i=1:size(reduced_data,1)
    entropy_reduced_data(i,:)=P'*reduced_data(i,:)';
end
syms beta1 beta2
beta=[beta1 beta2]';
% Class 0
x0=entropy_reduced_data(train_classess==0,:);
mu0=mean(x0)';
sigma0=cov(x0);
A0=-0.5*inv(sigma0);
b0=inv(sigma0)*mu0;
c0=-0.5*mu0'*sigma0*mu0-0.5*det(sigma0);
d0=beta'*A0*beta + b0'*beta + c0;

% Class 1
x1=entropy_reduced_data(train_classess==1,:);
mu1=mean(x1)';
sigma1=cov(x1);
A1=-0.5*inv(sigma1);
b1=inv(sigma1)*mu1;
c1=-0.5*mu1'*sigma1*mu1-0.5*det(sigma1);
d1=beta'*A1*beta + b1'*beta + c1;
boundary=d0-d1;
boundary_bayes_entropy=boundary;
boundary_euclidean_entropy=sqrt(sum(beta-mu0).^2)-sqrt(sum(beta-mu1).^2);
boundary_mah_entropy=-0.5*(beta-mu0)' * inv(sigma0) * (beta-mu0) + 0.5*(beta-mu1)' * inv(sigma1) * (beta-mu1);
figure
scatter(entropy_reduced_data(train_classess==0,1),entropy_reduced_data(train_classess==0,2))
hold on
scatter(entropy_reduced_data(train_classess==1,1),entropy_reduced_data(train_classess==1,2))
hold on
fimplicit(boundary==0)
legend('Class 1','Class 2','Boundary')
title('Feature Reduction by Using Minimum Entropy')

%% E: PCA
sigma=cov(reduced_data);
[V,lambda]=eig(sigma);
P=[V(:,3) V(:,2)];
P_pca=P;
for i=1:size(reduced_data,1)
    pca_reduced_data(i,:)=P'*reduced_data(i,:)';
end
syms beta1 beta2
beta=[beta1 beta2]';
% Class 0
x0=pca_reduced_data(train_classess==0,:);
mu0=mean(x0)';
sigma0=cov(x0);
A0=-0.5*inv(sigma0);
b0=inv(sigma0)*mu0;
c0=-0.5*mu0'*sigma0*mu0-0.5*det(sigma0);
d0=beta'*A0*beta + b0'*beta + c0;

% Class 1
x1=pca_reduced_data(train_classess==1,:);
mu1=mean(x1)';
sigma1=cov(x1);
A1=-0.5*inv(sigma1);
b1=inv(sigma1)*mu1;
c1=-0.5*mu1'*sigma1*mu1-0.5*det(sigma1);
d1=beta'*A1*beta + b1'*beta + c1;
boundary=d0-d1;
boundary_bayes_pca=boundary;
boundary_euclidean_pca=sqrt(sum(beta-mu0).^2)-sqrt(sum(beta-mu1).^2);
boundary_mah_pca=-0.5*(beta-mu0)' * inv(sigma0) * (beta-mu0) + 0.5*(beta-mu1)' * inv(sigma1) * (beta-mu1);
figure
scatter(pca_reduced_data(train_classess==0,1),pca_reduced_data(train_classess==0,2))
hold on
scatter(pca_reduced_data(train_classess==1,1),pca_reduced_data(train_classess==1,2))
hold on
fimplicit(boundary==0)
legend('Class 1','Class 2','Boundary')
title('Feature Reduction by Using PCA')
%% F: FLD
% Finding the best feature
for i=1:3
    selected_feature=reduced_data(:,i);
    p_w0=length(find(train_classess==0));
    p_w1=length(find(train_classess==1));
    
    mu_0=mean(selected_feature(train_classess==0,:));
    mu_0=mu_0';
    mu_1=mean(selected_feature(train_classess==1,:));
    mu_1=mu_1';
    
    mu_total=p_w0*mu_0+p_w1*mu_1;
    %mu_total=mu_0+mu_1;
    mu_total=mu_total';
    Sb=p_w0*(mu_0-mu_total)*(mu_0-mu_total)' + (mu_1-mu_total)*(mu_1-mu_total)';
    x0=selected_feature(train_classess==0,:);
    x1=selected_feature(train_classess==1,:);
    cov0=0;
    for j=1:length(x0)
        %cov0=cov0+(x0(j,:)'-mu_0)*(x0(j,:)'-mu_0)';
        cov0=cov0+(x0(j,:)'-mu_total)*(x0(j,:)'-mu_total)';
    end
    cov1=0;
    for j=1:length(x1)
        %cov1=cov1+(x1(j,:)'-mu_1)*(x1(j,:)'-mu_1)';
        cov1=cov1+(x1(j,:)'-mu_total)*(x1(j,:)'-mu_total)';
    end

    %Sw=p_w0*cov0+p_w1*cov1;
    Sw=cov0+cov1;
    J(i)=trace(Sb)/trace(Sw);
end
best_feature=find(J==max(J));

% Finding the second best feature
k=1:3;
k=setdiff(k,best_feature);
J=0;
for i=k
    selected_feature=[reduced_data(:,best_feature) reduced_data(:,i) ];
    p_w0=length(find(train_classess==0));
    p_w1=length(find(train_classess==1));
    
    mu_0=mean(selected_feature(train_classess==0,:));
    mu_0=mu_0';
    mu_1=mean(selected_feature(train_classess==1,:));
    mu_1=mu_1';
    
    mu_total=p_w0*mu_0+p_w1*mu_1;
    %mu_total=mu_0+mu_1;
    mu_total=mu_total';
    Sb=p_w0*(mu_0-mu_total)*(mu_0-mu_total)' + (mu_1-mu_total)*(mu_1-mu_total)';
    x0=selected_feature(train_classess==0,:);
    x1=selected_feature(train_classess==1,:);
    cov0=0;
    for j=1:length(x0)
        %cov0=cov0+(x0(j,:)'-mu_0)*(x0(j,:)'-mu_0)';
        cov0=cov0+(x0(j,:)'-mu_total)*(x0(j,:)'-mu_total)';
    end
    cov1=0;
    for j=1:length(x1)
        %cov1=cov1+(x1(j,:)'-mu_1)*(x1(j,:)'-mu_1)';
        cov1=cov1+(x1(j,:)'-mu_total)*(x1(j,:)'-mu_total)';
    end

    %Sw=p_w0*cov0+p_w1*cov1;
    Sw=cov0+cov1;
    J(i)=trace(Sb)/trace(Sw);
end
best_feature2=find(J==max(J));
fld_reduced_data=[reduced_data(:,best_feature) reduced_data(:,best_feature2)];
fld_best_featrues=[best_feature best_feature2];

syms beta1 beta2
beta=[beta1 beta2]';
% Class 0
x0=fld_reduced_data(train_classess==0,:);
mu0=mean(x0)';
sigma0=cov(x0);
A0=-0.5*inv(sigma0);
b0=inv(sigma0)*mu0;
c0=-0.5*mu0'*sigma0*mu0-0.5*det(sigma0);
d0=beta'*A0*beta + b0'*beta + c0;

% Class 1
x1=fld_reduced_data(train_classess==1,:);
mu1=mean(x1)';
sigma1=cov(x1);
A1=-0.5*inv(sigma1);
b1=inv(sigma1)*mu1;
c1=-0.5*mu1'*sigma1*mu1-0.5*det(sigma1);
d1=beta'*A1*beta + b1'*beta + c1;
boundary=d0-d1;
boundary_bayes_fld=boundary;
boundary_euclidean_fld=sqrt(sum(beta-mu0).^2)-sqrt(sum(beta-mu1).^2);
boundary_mah_fld=-0.5*(beta-mu0)' * inv(sigma0) * (beta-mu0) + 0.5*(beta-mu1)' * inv(sigma1) * (beta-mu1);
figure
scatter(fld_reduced_data(train_classess==0,1),fld_reduced_data(train_classess==0,2))
hold on
scatter(fld_reduced_data(train_classess==1,1),fld_reduced_data(train_classess==1,2))
hold on
fimplicit(boundary==0)
legend('Class 1','Class 2','Boundary')
title('Feature Reduction by Using FLD')

%% G
% Calculating Test Features
% Calculating the first feature: Band Power
band_f0=[0 4 8 14];
band_f1=[3 7 13 20];
for i=1:length(test_data)
    EEG=test_data{i};
    for j=1:6
        for k=1:4
            BW_Power_test(i,j,k)=bandpower(EEG(j,:),fs,[band_f0(k) band_f1(k)]);
        end
    end
end
% Calculating the second feature: Asymmetry Ratio
for i=1:length(test_data)
    c=1;
    for j=1:2:5
        for k=2:2:6
            for l=1:4
                R=BW_Power_test(i,k,l);
                L=BW_Power_test(i,j,l);
                AR_test(i,c)=(R-L)/(R+L);
                c=c+1;
            end
        end
    end
end
% Converting BW Power from rank 3 Tensor to rank 2 Tensor
for i=1:length(test_data)
    c=1;
    for j=1:6
        for k=1:4
            BW_Power_Feature_test(i,c)=BW_Power_test(i,j,k);
            c=c+1;
        end
    end
end

% Combining all features
all_test_features=[BW_Power_Feature_test AR_test];

% Extracting best features
reduced_test_features=all_test_features(:,all_best_features);

% Bayesian
% Bayesian Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_bayes_bayesian,[beta1 beta2 beta3],reduced_test_features(i,:));
    if dec>=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_bayes_bayesian=(TP+TN)/(TP+TN+FP+FN);
% Euclidean Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_euclidean_bayesian,[beta1 beta2 beta3],reduced_test_features(i,:));
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_euclidean_bayesian=(TP+TN)/(TP+TN+FP+FN);
% MAH Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_mah_bayesian,[beta1 beta2 beta3],reduced_test_features(i,:));
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_mah_bayesian=(TP+TN)/(TP+TN+FP+FN);

% Entropy
reduced_test_entropy=P_entropy'*reduced_test_features';
% Bayesian Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_bayes_entropy,[beta1 beta2],reduced_test_entropy(:,i)');
    if dec>=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_bayes_entropy=(TP+TN)/(TP+TN+FP+FN);
% Euclidean Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_euclidean_entropy,[beta1 beta2],reduced_test_entropy(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_euclidean_entropy=(TP+TN)/(TP+TN+FP+FN);
% MAH Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_mah_entropy,[beta1 beta2],reduced_test_entropy(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_mah_entropy=(TP+TN)/(TP+TN+FP+FN);

% PCA
reduced_test_pca=P_pca'*reduced_test_features';
% Bayesian Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_bayes_pca,[beta1 beta2],reduced_test_pca(:,i)');
    if dec>=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_bayes_pca=(TP+TN)/(TP+TN+FP+FN);
% Euclidean Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_euclidean_pca,[beta1 beta2],reduced_test_pca(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_euclidean_pca=(TP+TN)/(TP+TN+FP+FN);
% MAH Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_mah_pca,[beta1 beta2],reduced_test_pca(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_mah_pca=(TP+TN)/(TP+TN+FP+FN);
% FLD
reduced_test_fld=reduced_test_features(:,fld_best_featrues);
% Bayesian Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_bayes_fld,[beta1 beta2],reduced_test_fld(:,i)');
    if dec>=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_bayes_fld=(TP+TN)/(TP+TN+FP+FN);
% Euclidean Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_euclidean_fld,[beta1 beta2],reduced_test_fld(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_euclidean_fld=(TP+TN)/(TP+TN+FP+FN);
% MAH Distance
for i=1:size(reduced_test_features,1)
    dec=subs(boundary_mah_fld,[beta1 beta2],reduced_test_fld(:,i)');
    if dec<=0
        predicted_class(i)=0;
    else
        predicted_class(i)=1;
    end
end
TP=length(find(predicted_class==1 & test_classess==1));
TN=length(find(predicted_class==0 & test_classess==0));
FP=length(find(predicted_class==1 & test_classess==0));
FN=length(find(predicted_class==0 & test_classess==1));
ACC_mah_fld=(TP+TN)/(TP+TN+FP+FN);
















