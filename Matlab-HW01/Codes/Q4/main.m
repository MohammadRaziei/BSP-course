disp('Question 4 is running ...');
%% global variables or settings
addpath(fullfile('..','OSET'), fullfile('../OSET','Tools'), fullfile('../../../','DataSets'));%, fullfile('../OSET','Databases'))
addpath(fullfile('../OSET/Tools/Mex Files','Filters'),fullfile('../OSET/Tools/Mex Files','MedianFilter3'),fullfile('../OSET/Tools/Mex Files','TrimmedFilter3'))
addpath(fullfile('..','Functions'))

%% A
load('Signal1.mat');
fs = Signal1.fs;
ECG1 = Signal1.ECG - mean(Signal1.ECG);
figure; %PSD has been deprecated. Use PERIODOGRAM or PWELCH instead.
maskPlot = 1:2e3;
subplot(4,2,1);plot(Signal1.ECG(maskPlot));
subplot(4,2,2);pwelch(Signal1.ECG);

subplot(4,2,3);plot(ECG1(maskPlot));
subplot(4,2,4);pwelch(ECG1);
% b = TrimmedFilter(ECG1,'median',(.2*fs));
b=TrimmedFilter(ECG1,length(ECG1),'median',(.2*fs));
% b = TrimmedFilter(b,'median',(.4*fs));
b = TrimmedFilter(b,length(b),'median',(.4*fs));
b = LPFilter(b',10/fs);
ECG2 = ECG1 - b;
subplot(4,2,5);plot(ECG2(maskPlot));
subplot(4,2,6);pwelch(ECG2);
% Remove civil noise
[num,den] = iirnotch(60/(fs/2),0.25/(fs/2));
ECG3 = filter(num,den,ECG2);
subplot(4,2,7);plot(ECG3(maskPlot));
subplot(4,2,8);pwelch(ECG3);

%% B
peaksIdx = myFindPeaks2(ECG3);
RR = diff(peaksIdx) / fs;
% sdk = zeros(size(ECG3));
% sdk(peaksIdx) = ECG3(peaksIdx);
figure;
subplot(211);
hold on
plot(ECG3,'r');
% stem(sdk(maskPlot))
% plot(peaksIdx,RR,'*')

disp(all(peaksIdx' == myFindPeaks3(ECG3)));

idx2 = find(PeakDetection2(ECG3,180,1));
plot(idx2, ECG3(idx2),'*')
title('ECG Peak Detextion');
% pd = PeakDetection(ECG3,180,1);
% unique(pd)
subplot(212);
plot(RR);
hold on;
plot(RR,'o');
title('RR');
ax = gca;
ax.XTick = 1:23;
