function idx = myFindPeaks3(signalECG)

mainLobe = 3;
gaurd = 5;
Th = 100;

lag = floor(mainLobe/2)+gaurd;
lag2 = floor(lag/2)+1;
filt = [-1/(2*gaurd)*ones(1,gaurd) 1/mainLobe*ones(1,mainLobe) -1/(2*gaurd)*ones(1,gaurd)];

y = filter(filt,1,signalECG);
y = find(y > Th)-lag2;
diffY = diff(y);
NonOneY = find(diffY~=1);
idx = zeros(length(NonOneY)+1,1);
start = 1;
for i = 1:length(idx)-1
    endd = NonOneY(i);
    [~,v] = max(signalECG(y(start:endd)));
    idx(i) = y(start+v-1);
    start = endd+1;
end
[~,v] = max(signalECG(y(start:end)));
idx(i+1) = y(start+v-1);

end