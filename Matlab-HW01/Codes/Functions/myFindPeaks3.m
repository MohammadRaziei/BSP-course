function idx = myFindPeaks3(signalECG)

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
idx = zeros(length(NonOneY)+1,1);
start = 1;
for i = 1:length(idx)-1
    endd = NonOneY(i);
    start = y(start)-extra;
    endd = y(endd);
    
    [~,v] = max(signalECG(start:endd));
    idx(i) = (start+v-1);
    start = NonOneY(i)+1;
end
% [~,v] = max(signalECG(start:y(end)));
% idxx(i+1) = start+v-1;
% start = y(end)-extra;
start = y(start)-extra;
[~,v] = max(signalECG(start:y(end)));
idx(i+1) = start+v-1;



end