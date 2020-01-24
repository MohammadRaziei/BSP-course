disp('Question 1 is running ...');
%% global variables or settings
addpath(fullfile('../../../','DataSets'));
load('eegdata.mat')
fs = 250;
t = ((1:2500)-1)/fs;
task = data{1}{4};
plot(t,task(1,:))
title(['C3',',',data{1}{1},',',data{1}{2},',',data{1}{3}]);xlabel('t(sec)');

C3=1;C4=2;P3=3;P4=4;O1=5;O2=6;EOG=7;

% figure('units','normalized','outerposition',[0 0 1 1])
% annotation('line', 0.37*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);
% annotation('line', 0.65*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);

%% A


