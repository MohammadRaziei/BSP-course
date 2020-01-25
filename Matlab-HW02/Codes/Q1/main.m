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
Noise = 1*cos(2*pi*60*t + 2*pi*rand);

EEG = double(task(c3,:));
% EEG_N = EEG + Noise;



% alf = dsp.AdaptiveLatticeFilter;
% [y,err] = alf(x,d);
% dsp.RLSFilter

% n0 = 110;
% d = EEG_N(1:end-n0);
% x = EEG_N(n0+1:end);

d = EEG;
rlsFilt = dsp.RLSFilter(10);
[y,e] = rlsFilt(Noise',d');
% lmsFilt = dsp.LMSFilter;
% [y,e] = lmsFilt(x',d');
EEG_DN = d-y';

%% B

figure;
m = 1;
subplot(211);plot(EEG(m:end),'k');
title('Raw EEG');
% hold on; plot(EEG(m:end));
subplot(212);plot(EEG_DN(m:end));
title('Denoised EEG');
ylim(20*[-1,1])



figure;
subplot(211); pwelch(EEG,fs);
subplot(212); pwelch(EEG_DN,fs);

%% C
SNR_in= @(s,n) 10*log10(sqrt(sum(s.^2))/sqrt(sum(n.^2)));
SNR_out= @(s,s_hat) 10*log10(sqrt(sum(s.^2))/sqrt(sum((s-s_hat).^2)));
SNR_imp = @(s,n,s_hat) SNR_out(s,s_hat) - SNR_in(s,n);


s=d;
s_hat=d-y';
snr_in = SNR_in(s, Noise)
snr_out = SNR_out(s, s_hat)
snr_imp = SNR_imp(s,Noise,s_hat)

%% D
A = 0.5:0.5:5;
snrInArr = zeros(1,length(A));
snrOutArr = zeros(1,length(A));
snrImpArr = zeros(1,length(A));
for i=1:length(A)
    n = A(i)*cos(2*pi*60.*t + (2*pi)*rand);
    lmsFilt=dsp.LMSFilter(10);
    lmsFilt.StepSize=0.003;
    [y,e,wts] = lmsFilt(n',d');
    
    s=d;
    s_hat=d-y';
    
    snrInArr(i) = SNR_in(s,n);
    snrOutArr(i) = SNR_out(s,s_hat);
    snrImpArr(i)= SNR_imp(s,n,s_hat);
end
figure
plot(snrInArr,snrImpArr)
%% E


w0=2*pi*60/fs;















