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
signalECG = ECG;
% maxTh = max(ECG);

% buffLen = 120;
% signalECGtemp = reshape(signalECG,[],buffLen);
% len = length(signalECGtemp);
% [v,idx] = max(signalECGtemp);
% 
% idx = idx + (0:buffLen-1)*len;
% idx(find(signalECG(idx) < 0.65*maxTH)) = [];

mainLobe = 3;
gaurd = 5;
th = 75;
extra = 8;

% lag = floor(mainLobe/2)+gaurd;
% lag2 = ceil(lag/2)+1;
filt = [-1/(2*gaurd)*ones(1,gaurd) 1/mainLobe*ones(1,mainLobe) -1/(2*gaurd)*ones(1,gaurd)];

y0 = filter(filt,1,signalECG);
y = find(y0 > th);
diffY = diff(y);
NonOneY = find(diffY~=1);
idxx = zeros(length(NonOneY)+1,1);
start = 1;
for i = 1:length(idxx)-1
    endd = NonOneY(i);
    start = y(start)-extra;
    endd = y(endd);
    
    [~,v] = max(signalECG(start:endd));
    idxx(i) = (start+v-1);
    start = NonOneY(i)+1;
end
% [~,v] = max(signalECG(start:y(end)));
% idxx(i+1) = start+v-1;
% start = y(end)-extra;
start = y(start)-extra;
[~,v] = max(signalECG(start:y(end)));
idxx(i+1) = start+v-1;

% 
clf;
plot(idxx,signalECG(idxx),'*');
% plot(y)
hold on;
plot(y0,'y')
plot(signalECG)
