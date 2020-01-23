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
subplot(121);
plot(Error)
title("levinson")
subplot(122);
plot(AIC)
title("AIC")
% [a,e] = levinson(r,n);

%% D
[Pxx,w_D] = pyulear(Y,5);
figure;hold on
plot(w(1:1000)/pi,20*log10(abs(h(1:1000))))
plot(w_D/pi,10*log10(Pxx))
legend('True Power Spectral Density','pyulear PSD Estimate')
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% E

% approach — Algorithm for computing AR model
% 'fb' (default) | 'burg' | 'gl' | 'ls' | 'yw'
% Algorithm for computing the AR model, specified as one of the following values:
% 
% 'burg': Burg's lattice-based method. Solves the lattice filter equations using the harmonic mean of forward and backward squared prediction errors.
% 'fb': (Default) Forward-backward approach. Minimizes the sum of a least-squares criterion for a forward model, and the analogous criterion for a time-reversed model.
% 'gl': Geometric lattice approach. Similar to Burg's method, but uses the geometric mean instead of the harmonic mean during minimization.
% 'ls': Least-squares approach. Minimizes the standard sum of squared forward-prediction errors.
% 'yw': Yule-Walker approach. Solves the Yule-Walker equations, formed from sample covariances.

ar_alg = ["burg","fb","gl","ls","yw"];
ar_A = [];
for i = 1:length(ar_alg)
    sys = ar(data,5,ar_alg(i));
    ar_A = [ar_A;sys.A];
end
%%
figure;hold on    
plot(w(1:1e3)/pi,20*log10(abs(h(1:1e3))),'r');
legend_cell = cell(1,length(ar_alg)+2);
legend_cell{1} = 'True Power Spectral Density';
for i = 1:length(ar_alg)
    [h_E,w_E] = freqz(1,ar_A(i,:),'whole',2001);
    plot(w_E(1:1e3)/pi,20*log10(abs(h_E(1:1e3))));
    legend_cell{i+1} = ['AR Power Estimate with "' char(ar_alg(i)) '"'];
end
plot(w_D/pi,10*log10(Pxx))
legend_cell{i+2} = 'pyulear PSD Estimate';
legend(legend_cell)
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
ylim([-50,100])
