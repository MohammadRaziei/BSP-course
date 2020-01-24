disp('Question 5 is running ...');close all;
%% global variables or settings
% addpath(fullfile('..','OSET'), fullfile('../OSET','Tools');
% addpath(fullfile('../OSET/Tools/Mex Files','Filters'),fullfile('../OSET/Tools/Mex Files','MedianFilter3'),fullfile('../OSET/Tools/Mex Files','TrimmedFilter3'))
addpath(fullfile('../../../','DataSets'));
addpath(fullfile('..','Functions'))

%% A
load('18184m.mat');

% plot(t,ECG);

% resample(ones(1,5),3,5)
ECG = val(1,:);
% Data of 30 seconds
ECG=ECG(1,1:length(ECG)/2);
fs=128;
t=0:1/fs:(length(ECG)-1)/fs;

%% B
% peaksIdx =find(PeakDetection2(ECG,100));
peaksIdx = myFindPeaks3(ECG);

figure;
plot(t,ECG);hold on;
plot(peaksIdx/fs,ECG(peaksIdx),'*')
title('ECG Signal')
xlabel('Time (s)')
%% C
SNR=15:-5:5;
ECG_N = zeros(length(SNR),length(ECG));

figure;
subplot(length(SNR)+1,1,1);plot(t,ECG,'k');
title('Primery ECG')
for i=1:length(SNR)
    % Adding white guassian noise
    ECG_N(i,:)=awgn(ECG,SNR(i),'measured');
    subplot(length(SNR)+1,1,i+1);plot(t,ECG_N(i,:));
    title(['Noisy ECG with SNR=' num2str(SNR(i))])
end
xlabel('Time (s)')


%% D

numSegments=length(peaksIdx)-1;

lag = floor(0.6*min(diff(peaksIdx)));

figure;
for channel = 1:length(SNR)
    % Algorithm 1:
    beatLen = floor(mean(diff(peaksIdx)));
    SynchSignal=zeros(numSegments,beatLen);
    for i=1:numSegments-1
        RR = peaksIdx(i)-lag:peaksIdx(i+1)-lag;
        SynchSignal(i,:)=resample(ECG_N(channel,RR),beatLen,length(RR));
    end
    SynchAvg = mean(SynchSignal);
    subplot(length(SNR),1,channel) ;plot(SynchAvg)
    title(['SNR = ' num2str(SNR(channel))])
end
figure;
for channel = 1:length(SNR)
    % Algorithm 2:
    beatLen = max(diff(peaksIdx))+1;
    % lag = floor(0.5*min(diff(peaksIdx)));
    SynchAvg=zeros(1,beatLen);
    for i=1:numSegments
        RR = peaksIdx(i)-lag:peaksIdx(i+1)-lag;
        temp = (1/numSegments)*ECG_N(channel,RR);
        SynchAvg(1:length(RR)) = SynchAvg(1:length(RR)) + temp;
    end
    subplot(length(SNR),1,channel) ;plot(SynchAvg)
    title(['SNR = ' num2str(SNR(channel))])
end




