disp('Question 4 is running ...');
%% global variables or settings
addpath(fullfile('..','OSET'), fullfile('../OSET','Tools'), fullfile('../../../','DataSets'));%, fullfile('../OSET','Databases'))
addpath(fullfile('../OSET/Tools/Mex Files','Filters'),fullfile('../OSET/Tools/Mex Files','MedianFilter3'),fullfile('../OSET/Tools/Mex Files','TrimmedFilter3'))
%% A
load('Signal1.mat');
fs = Signal1.fs;
ECG1 = Signal1.ECG - mean(Signal1.ECG);
figure;
subplot(121);plot(ECG1);
subplot(122);pwelch(ECG1);
% b = TrimmedFilter(ECG1,'median',(.2*fs));
% b = TrimmedFilter(b,'median',(.4*fs));
% b = LPFilter(b,10/fs);
b = mean(ECG1);
ECG2 = ECG1 - b;

% Remove civil noise
[num,den] = iirnotch(60/(fs/2),0.25/(fs/2));
ECG3 = filter(num,den,ECG2);
figure; pwelch(ECG3);