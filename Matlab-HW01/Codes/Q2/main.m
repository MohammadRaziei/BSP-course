disp('Question 2 is running ...');
%% global variables or settings
A = [1 -1.1 -0.0975 0.019 1.0825 -0.904];
B = [1 0.5];

sys = tf(B,A,1);
% sys = zpk(B,A,1);
pzmap(sys)
% model = arima('Constant',0,'AR',{1.1,0.0975,-0.019,-1.0825,0.904},'MA',{0.5},'Variance',1);
%% A
[h,w] = freqz(B,A,'whole',2001);
figure
plot(w/pi,20*log10(abs(h)))
ax = gca;
% ax.YLim = [-100 20];
ax.XTick = 0:.2:2;
% ax.XTick = [0 0.2148 0.75 1.7852 1.25 2];

xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
% armax(data,

%% B
Num = 10000;
U = wgn(1,Num,1);
% hist(U(1,:))
Y = filter(B,A, U);
figure;
plot(Y);

%% C
% close all; 
data = iddata(Y');
% sys0 = idpoly(A,B);

[r,lg] = xcorr(Y,'biased');
r(lg<0) = [];

Error = []; AIC = [];
for i = 1:10
    [a,e] = levinson(r,i);
    Error = [Error, e];
    sys = ar(data,i);
    aic_raw = aic(sys,'aic');
    AIC = [AIC, aic_raw];
end
figure;
plot(Error)
title("levinson")
figure;
plot(AIC)
title("AIC")
% [a,e] = levinson(r,n);

%% D
% sys = ar(data,5);
[Pxx,w_D] = pyulear(Y,5);
figure;hold on
plot(w(1:1000)/pi,20*log10(abs(h(1:1000))))
plot(w_D/pi,10*log10(Pxx))
legend('True Power Spectral Density','pyulear PSD Estimate')

%% E


%% F