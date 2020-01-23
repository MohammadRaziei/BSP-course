function idx = myFindPeaks2(signalECG)

maxTh = max(signalECG);

buffLen = 120;
signalECGtemp = reshape(signalECG,[],buffLen);
len = length(signalECGtemp);
[~,idx] = max(signalECGtemp);

idx = idx + (0:buffLen-1)*len;
idx(find(signalECG(idx) < 0.65*maxTh)) = [];

end