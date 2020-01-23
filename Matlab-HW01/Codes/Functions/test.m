% diffSig = diff(ECG3);
% thy = 80;
% thx = 300;
% peaks1 = find(diffSig > thy);
% diffpeaks = diff(peaks1);
% peaks1(find(diffpeaks < thx)+1) = [];
% 
% peaks2 = find(diffSig < -thy);
% diffpeaks1 = -diff(peaks2(end:-1:1));
% peaks2(length(peaks2) - find(diffpeaks1 < thx)) = [];
% 
% % peaksIdx = floor(0.5*(peaks1+peaks2));
signalECG = ECG3;
maxTh = max(ECG3);

% buffLen = 120;
% signalECGtemp = reshape(signalECG,[],buffLen);
% len = length(signalECGtemp);
% [v,idx] = max(signalECGtemp);
% 
% idx = idx + (0:buffLen-1)*len;
% idx(find(signalECG(idx) < 0.65*maxTH)) = [];

mainLobe = 3;
gaurd = 5;
lag = floor(mainLobe/2)+gaurd;
lag2 = floor(lag/2)+1;
filt = [-1/(2*gaurd)*ones(1,gaurd) 1/mainLobe*ones(1,mainLobe) -1/(2*gaurd)*ones(1,gaurd)];

y = filter(filt,1,signalECG);
y = find(y > 100)-lag2;
diffY = diff(y);
NonOneY = find(diffY~=1);
idxx = zeros(length(NonOneY)+1,1);
start = 1;
for i = 1:length(idxx)-1
    endd = NonOneY(i);
    [~,v] = max(signalECG(y(start:endd)));
    idxx(i) = y(start+v-1);
    start = endd+1;
end
[~,v] = max(signalECG(y(start:end)));
idxx(i+1) = y(start+v-1);

% 
% clf;
% plot(y,ECG3(y),'*');
% % plot(y)
% hold on;
% plot(ECG3)