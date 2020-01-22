disp('Question 3 is running ...');
%% global variables or settings
A = [1 -1.1 -0.0975 0.019 1.0825 -0.904];
B = [1 0.5]; 

%% A
Num = 1000;
U = wgn(1,Num,1);
Y = filter(B,A,U);

[Ry_b,~] = xcorr(Y,'biased');
[Ry_u,~] = xcorr(Y,'unbiased');
% r(lg<0) = [];

Sy_BT_b = fft(Ry_b);
Sy_BT_u = fft(Ry_u);
[h,w] = freqz(B,A,'whole',2*Num-1);
% w = linspace(0,1,2*Num-1);
figure; 
subplot(211); hold on;
plot(w/pi, 10*log10(abs(Sy_BT_b)));
plot(w/pi, 10*log10(abs(Sy_BT_u)));
plot(w/pi, 20*log10(abs(h)));
ax = gca;
ax.XTick = 0:.2:2;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
legend("biased","unbiased","freqz");


%% B
Sy_BT_b = zeros(1,2*Num-1);
Sy_BT_u = zeros(1,2*Num-1);
for i = 1:100
    Yu = filter(B,A,wgn(1,Num,1));
    [Ry_b,~] = xcorr(Yu,'biased');
    [Ry_u,~] = xcorr(Yu,'unbiased');
    Sy_BT_b = Sy_BT_b + fft(Ry_b)/100;
    Sy_BT_u = Sy_BT_u + fft(Ry_u)/100;
end
subplot(212); hold on;
plot(w/pi, 10*log10(abs(Sy_BT_b)));
plot(w/pi, 10*log10(abs(Sy_BT_u)));
plot(w/pi, 20*log10(abs(h)));
ax = gca;
ax.XTick = 0:.2:2;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
legend("biased","unbiased","freqz");

%% C
win = hanning(length(Y));
[pxx,w0] = periodogram(Y,win,length(Y));%,'power');%,'reassigned');
% [pxx,w0] = periodogram(Y,'power');
% [pxx,w0] = periodogram(Y);

figure; 
hold on;
plot(w0/pi, 10*log10(pxx));
plot(w(1:1000)/pi, 20*log10(abs(h(1:1000))));
ax = gca;
ax.XTick = 0:.2:2;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
legend("periodogram","freqz");

%% D
figure; 
p = [0.25,0.50,0.75];
winn = [];
title_winn = 'rectangular window';
for ii = 1:2
    subplot(1,2,ii); hold on;  
    if ii == 2
        winn = win;
        title_winn = 'hanning window';
    end
    legend_cell = cell(1,length(p)+1);
    for i = 1:length(p)
        [pyy,w1] = pwelch(Y,winn,p(i)*length(Y));
        plot(w1/pi, 10*log10(pyy)); 
        legend_cell{i} = ['pwelch for ' num2str(p(i)*100) '%'];
    end
    plot(w(1:1000)/pi, 20*log10(abs(h(1:1000))));
    legend_cell{i+1} = 'freqz';
    ax = gca;
    ax.XTick = 0:.2:2;
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude (dB)');
    title(title_winn);
    legend(legend_cell);
end
%% E
data = iddata(Y');
% econometricModeler
sys = armax(data,[5 1]);

[H_arma,w_arma] = freqz(sys.C,sys.A,'whole',2*Num-1);

figure
hold on;
plot(w_arma/pi,20*log10(abs(H_arma)));
plot(w/pi,20*log10(abs(h)));
ax = gca;
% ax.YLim = [-100 20];
ax.XTick = 0:.2:2;
% ax.XTick = [0 0.2148 0.75 1.7852 1.25 2];

xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
legend("ARMA(5,1)","Real spectrum");

