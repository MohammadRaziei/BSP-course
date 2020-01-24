disp('Question 4 is running ...');
%% global variables or settings
addpath(fullfile('..','OSET'), fullfile('../OSET','Tools'));
addpath(fullfile('../OSET/Tools/Mex Files','Filters'),fullfile('../OSET/Tools/Mex Files','MedianFilter3'),fullfile('../OSET/Tools/Mex Files','TrimmedFilter3'))
addpath(fullfile('../../../','DataSets'));
addpath(fullfile('..','Functions'))

%% A
load('Signal1.mat');
ECG1 = Signal1.ECG - mean(Signal1.ECG);
fs = Signal1.fs;
t = (0:length(ECG1)-1)/fs;
figure('units','normalized','outerposition',[0 0 1 1]);
%PSD has been deprecated. Use PERIODOGRAM or PWELCH instead.
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
figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2,2,[1,2]);
subplot(211)
hold on
plot(t,ECG3,'r');
% stem(sdk(maskPlot))
% plot(peaksIdx,RR,'*')

disp(all(peaksIdx' == myFindPeaks2(ECG3)));

idx = find(PeakDetection(ECG3,1/fs));
idx2 = find(PeakDetection2(ECG3,180,1));

plot(idx2/fs,ECG3(idx2),'*')
title('ECG Peak Detextion');


% subplot(223);
subplot(212);
% plot(1+(0:length(RR_HR)-1)/5,RR_HR);
plot(RR_idx,RR);
hold on;
plot(RR_idx,RR,'o');
title('RR');
% ax = gca;
% ax.XTick = 1:2:length(RR);
xlabel('indices');

% subplot(2,2,4);
% [w,pRR] = 
% pwelch(RR);
% RR = RR_HR;
%% C
[b,a]=butter(5,6/(fs/2),'high');
RR_filtered=filtfilt(b,a,RR);

figure('units','normalized','outerposition',[0 0 1 1]); 

% Periodogram
pRR = abs(fft(RR_filtered)).^2;
w = (0:length(RR)-1)/(length(RR)/2);
mask = 1:length(w)/2;
subplot(311);plot(w(mask),pRR(mask))
title('Periodogram')
% BT
pRR = abs(fft(xcorr(RR_filtered,'biased'))).^2;
w = (0:length(pRR)-1)/(length(pRR)/2);
mask = 1:length(w)/2;
subplot(312);plot(w(mask),pRR(mask))
title('BT')
% Welch
[pRR,w] = pwelch(RR_filtered);
subplot(313);plot(w/pi,pRR)
title('Welch')

% AR
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
p = [2:10 20];
% legend_cell = cell(1,length(p));
for i = 1:length(p)
    [pRR,w] = pyulear(RR_filtered,p(i));
    subplot(ceil(length(p)/2),2,i);plot(w/pi,abs(pRR))
%     legend_cell{i} = ['AR Power Estimate with p = ' num2str(p(i))];
    title(['AR Power Estimate with p = ' num2str(p(i))]);
end
% title('AR')
% legend(legend_cell)
% xlabel('Normalized Frequency (\times\pi rad/sample)')

% ylim([-50,100])

for i = 1:length(p)
    figure('units','normalized','outerposition',[0 0.2 0.8 0.8]);
    subplot(121);
    sys = ar(RR_filtered,p(i));
    [h,w]=freqz(1,sys.A);
    plot(w,abs(h).^2)
    title(['Predicted PSD by AR(',num2str(p(i)),')'])
    subplot(122);
    pzmap(sys)
    title(['Pole-Zero Map for AR(',num2str(p(i)),') Model'])
end


