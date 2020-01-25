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
% Noise = 3*cos(2*pi*60*t + 2*pi*rand(1,length(t)));

% EEG = double(task(c3,:));
% EEG_N = EEG + Noise;



% alf = dsp.AdaptiveLatticeFilter;
% [y,err] = alf(x,d);
% dsp.RLSFilter

n0 = 110;
% d = EEG_N(1:end-n0);
% x = EEG_N(n0+1:end);
d = EEG_N;
x = Noise;
rlsFilt = dsp.RLSFilter;
[y,e] = rlsFilt(x',d');
% lmsFilt = dsp.LMSFilter;
% [y,e] = lmsFilt(x',d');
EEG_DN = e';

%% B

figure;
m = 1;
subplot(311);plot(EEG(m:end),'k');
title('primary EEG');
subplot(312);plot(EEG_N(m:end),'g');
title('noisy EEG');
% hold on; plot(EEG(m:end));
subplot(313);plot(EEG_DN(m:end));
title('denoised noisy EEG');
ylim(20*[-1,1])



figure;
subplot(311); pwelch(EEG,fs);
subplot(312); pwelch(EEG_N,fs);
subplot(313); pwelch(EEG_DN,fs);

%% C
SNR_in= @(s,n) 10*log10(sqrt(sum(s.^2))/sqrt(sum(n.^2)));
SNR_out= @(s,y) 10*log10(sqrt(sum(s.^2))/sqrt(sum((s - y).^2)));
SNR_imp = @(SNRout,SNRin) SNRout-SNRin;


s=d;
s_hat=d-y';
snr_in = SNR_in(s, Noise)
snr_out = SNR_out(s, s_hat)
snr_imp = SNR_imp(snr_out,snr_in)

%% D
A=[0.5:0.5:5];
for i=1:length(A)
    phi=(2*pi)*rand;
    fs=250;
    dt=1/fs;
    t=0:dt:(2500-1)*dt;
    n=A(i)*cos(2*pi*60.*t + phi);

    L=10;
    rls=dsp.LMSFilter(L);
    rls.StepSize=0.003;
    [n_hat,err,wts] = rls(n',(d)');
    s=d;
    s_hat=d-n_hat';
    SNR_in=10*log10(sqrt(sum(s.^2))/sqrt(sum(n.^2)));
    SNR_in_vec(i)=SNR_in;
    
    SNR_out=10*log10(sqrt(sum(s.^2))/sqrt(sum((s-s_hat).^2)));
    SNR_imp=SNR_out-SNR_in;
    SNR_imp_vec(i)=SNR_imp;
end
figure
plot(SNR_in_vec,SNR_imp_vec)









