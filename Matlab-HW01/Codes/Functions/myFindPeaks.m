function peaksIdx = myFindPeaks(signalECG)

diffSig = diff(signalECG);
thy = 80;
thx = 300;
peaks1 = find(diffSig > thy);
diffpeaks = diff(peaks1);
peaks1(find(diffpeaks < thx)+1) = [];

peaks2 = find(diffSig < -thy);
diffpeaks1 = -diff(peaks2(end:-1:1));
peaks2(length(peaks2) - find(diffpeaks1 < thx)) = [];

peaksIdx = ceil(0.5*(peaks1+peaks2));

end