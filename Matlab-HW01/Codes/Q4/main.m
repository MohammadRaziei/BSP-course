disp('Question 4 is running ...');
%% global variables or settings
addpath(fullfile('..','OSET'), fullfile('../OSET','Tools'), fullfile('../../../','DataSets'));%, fullfile('../OSET','Databases'))
addpath(fullfile('../OSET/Tools/Mex Files','Filters'),fullfile('../OSET/Tools/Mex Files','MedianFilter3'),fullfile('../OSET/Tools/Mex Files','TrimmedFilter3'))
addpath(fullfile('..','Functions'))

%% A
load('Signal1.mat');
ECG1 = Signal1.ECG - mean(Signal1.ECG);
fs = Signal1.fs;
t = (0:length(ECG1)-1)/fs;
figure; %PSD has been deprecated. Use PERIODOGRAM or PWELCH instead.
maskPlot = 1:2e3;
subplot(4,2,1);plot(t(maskPlot),Signal1.ECG(maskPlot));title('ECG');
subplot(4,2,2);pwelch(Signal1.ECG);

subplot(4,2,3);plot(t(maskPlot),ECG1(maskPlot));title('ECG1');
subplot(4,2,4);pwelch(ECG1);
% b = TrimmedFilter(ECG1,'median',(.2*fs));
b=TrimmedFilter(ECG1,length(ECG1),'median',(.2*fs));
% b = TrimmedFilter(b,'median',(.4*fs));
b = TrimmedFilter(b,length(b),'median',(.4*fs));
b = LPFilter(b',10/fs)';
ECG2 = ECG1 - b;
subplot(4,2,5);plot(t(maskPlot),ECG2(maskPlot));title('ECG2');

subplot(4,2,6);pwelch(ECG2);
% Remove civil noise
[num,den] = iirnotch(60/(fs/2),0.25/(fs/2));
ECG3 = filter(num,den,ECG2);
subplot(4,2,7);plot(t(maskPlot),ECG3(maskPlot));title('ECG3');
subplot(4,2,8);pwelch(ECG3);

%% B
peaksIdx = myFindPeaks3(ECG3);
RR = diff(peaksIdx) / fs;
RR_HR = resample(RR,5,1);
RR_HR = RR_HR(1:end-4);
RR_idx = 1:length(RR);

% sdk = zeros(size(ECG3));
% sdk(peaksIdx) = ECG3(peaksIdx);
figure;
subplot(2,2,[1,2]);
hold on
plot(t,ECG3,'r');
% stem(sdk(maskPlot))
% plot(peaksIdx,RR,'*')

disp(all(peaksIdx' == myFindPeaks2(ECG3)));

idx2 = find(PeakDetection2(ECG3,180,1));

plot(idx2/fs,ECG3(idx2),'*')
title('ECG Peak Detextion');


% subplot(223);
subplot(211);
% plot(1+(0:length(RR_HR)-1)/5,RR_HR);
plot(RR_idx,RR);
hold on;
plot(RR_idx,RR,'o');
title('RR');
% ax = gca;
% ax.XTick = 1:2:length(RR);
xlabel('indices');

% subplot(2,2,4);
subplot(212);
% [w,pRR] = 
% pwelch(RR);
% RR = RR_HR;
pRR = fft(RR).^2;
w = (0:length(RR)-1)/length(RR);
plot(w,10*log10(abs(pRR)))

%% D




