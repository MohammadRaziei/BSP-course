disp('Question 1 is running ...');
%% global variables or settings
close all;
addpath(fullfile('../../../','DataSets'));
load('eegdata.mat')
fs = 250;
t = ((1:2500)-1)/fs;
task = data{1}{4};
plot(t,task(1,:))
title(['C3',',',data{1}{1},',',data{1}{2},',',data{1}{3}]);xlabel('t(sec)');

c3=1;c4=2;p3=3;p4=4;o1=5;o2=6;eog=7;

% figure('units','normalized','outerposition',[0 0 1 1])
% annotation('line', 0.37*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);
% annotation('line', 0.65*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);

%% A
Noise = 3*cos(2*pi*60*t + randn(1,length(t)));

EEG = task(c3,:);
EEG_N = EEG + Noise;



% alf = dsp.AdaptiveLatticeFilter;
% [y,err] = alf(x,d);
% dsp.RLSFilter

n0 = 110;
d = EEG_N(1:end-n0);
x = EEG_N(n0+1:end);
rlsFilt = dsp.RLSFilter;
[y,e] = rlsFilt(x,d);
% lmsFilt = dsp.LMSFilter;
% [y,e] = lmsFilt(x',d');

figure;
m = 1;
subplot(311);plot(EEG(m:end));
subplot(312);
plot(EEG_N(m:end),'g');
% hold on; plot(EEG(m:end));
subplot(313);plot(e(m:end));
ylim(20*[-1,1])



figure;
m = 1;
subplot(311); pwelch(EEG,fs);
subplot(312); pwelch(EEG_N,fs);
subplot(313); pwelch(e,fs);

