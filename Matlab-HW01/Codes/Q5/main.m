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
t = (0:length(ECG)-1)/fs;

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
    SynchSignal=zeros(numSegments-1,beatLen);
    for i=2:numSegments
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

% Robust Weighted Averaging
figure
SynchSignal=zeros(numSegments,beatLen);
for channel = 1:length(SNR)
    for i=1:numSegments
        RR = peaksIdx(i)-lag:peaksIdx(i+1)-lag; RRLen = length(RR);
        SynchSignal(i,1:RRLen) = ECG_N(channel,RR);
    end
    m=5; 
    % Initializing V
    v=zeros(1,beatLen);
    for k =1:100
        w = zeros(1,numSegments);
        for i=1:numSegments
            num=(norm(SynchSignal(i,:)-v))^(2/(1-m));
            den=0;
            for j=1:numSegments
                den = den+(norm(SynchSignal(j,:)-v))^(2/(1-m));
            end
            w(i)=num/den;
        end
        num=0;
        den=0;
        for i=1:numSegments
            num=num+(w(i)^m)*SynchSignal(i,:);
            den=den+(w(i)^m);
        end
        v=num/den;
    end
subplot(length(SNR),1,channel)
plot(v)
title(['Result of Robust Weighted Averaging for Signal with SNR=',num2str(SNR(channel))])
end

%% E
load('118e00m.mat');

NoisyECG  = val(2,:); % Data of 1 minute
NoisyECG = NoisyECG(1,1:length(NoisyECG)/2); % Data of 30 seconds
fs=360;
t = (0:length(NoisyECG)-1)/fs;

peaksIdx = myFindPeaks3(NoisyECG);

figure; hold on
plot(peaksIdx/fs, NoisyECG(peaksIdx),'*')

plot(t,NoisyECG)
title('ECG Signal')
xlabel('Time (s)')




% Synchronous Averaging
lag = floor(0.5*min(diff(peaksIdx)));
numSegments=length(peaksIdx) - 1;
beatLen = 2 + max(diff(peaksIdx));
SynchAvg=zeros(1,beatLen);
SynchAvg2=zeros(1,beatLen);
% lag = 0;
for i=2:numSegments
    RR= peaksIdx(i)-lag:peaksIdx(i+1)-lag; RRLen = length(RR);
    SynchAvg(1:RRLen)= SynchAvg(1:RRLen)+ (1/(numSegments-1))* NoisyECG(RR);
    SynchAvg2 = SynchAvg2 + (1/(numSegments-1))* resample(NoisyECG(RR),beatLen,RRLen);

end
% SynchAvg = fftshift(SynchAvg);
figure
subplot(211);plot(SynchAvg)
title('Result of Synchronous Averaging for Noise Stress Test Database - Algorithm 1')
subplot(212);plot(SynchAvg2)
title('Result of Synchronous Averaging for Noise Stress Test Database - Algorithm 2')

% Robust Weighted Averaging

SynchSignal=zeros(numSegments-2,beatLen);
for i=2:numSegments
    RR=peaksIdx(i)-lag:peaksIdx(i+1)-lag;
    SynchSignal(i-1,1:length(RR))= NoisyECG(RR);
end
v=zeros(1,beatLen);
for k=1:10
    w = zeros(1,numSegments-1);
    for i=1:numSegments-1
        num=(norm(SynchSignal(i,:)-v))^(2/(1-m));
        den=0;
        for j=1:numSegments-1
            den=den+(norm(SynchSignal(j,:)-v))^(2/(1-m));
        end
        w(i)=num/den;
    end
    num=0;
    den=0;
    for i=1:numSegments-1
        num=num+(w(i)^m)*SynchSignal(i,:);
        den=den+(w(i)^m);
    end
    v=num/den;
end
figure
plot(v)
title('Result of Robust Weighted Averaging for Noise Stress Test Database')


