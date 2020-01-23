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
% b = TrimmedFilter(b,'median',(.4*fs));
% b = LPFilter(b,10/fs);
b = mean(ECG1);
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
% sdk = zeros(size(ECG3));
% sdk(peaksIdx) = ECG3(peaksIdx);
figure;hold on
plot(ECG3,'r');
% stem(sdk(maskPlot))
plot(peaksIdx,ECG3(peaksIdx),'*')








