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
% signalTh = (ECG3>0.65*maxTh) + 0;
% th_min = min(find(diff(signalTh) == 1));
% th_mean = mean(diff(find(diff(signalTh) == 1)));
% th = floor(th_min*0.25 + th_mean*0.75);

buffLen = 120;
signalECGtemp = reshape(signalECG,[],buffLen);
len = length(signalECGtemp);
[v,idx] = max(signalECGtemp);

idx = idx + (0:buffLen-1)*len;
idx(find(signalECG(idx) < 0.65*maxTH)) = [];
