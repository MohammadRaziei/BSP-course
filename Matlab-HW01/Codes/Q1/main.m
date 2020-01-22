disp('Question 1 is running ...');
%% global variables or settings
Num = 8000;

figure('units','normalized','outerposition',[0 0 1 1])
annotation('line', 0.37*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);
annotation('line', 0.65*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);

%% A
Noise = 3*(rand(Num,1)-0.5);
subplot(331)
plot(Noise);
title("A")
ylim([-5,5]);
subplot(334)
hist(Noise)
subplot(337)
pwelch(Noise)

%% B
Noise = randn(Num,1);
subplot(332)
plot(Noise)
title("B")
ylim([-5,5]);
subplot(335)
hist(Noise)
subplot(338)
pwelch(Noise)

%% C
a = rand(Num,1); b = rand(Num,1);
randNoise = @(var,mu) var*cos(2*pi*b).*sqrt(-2*log(1-a))+mu;
Noise = randNoise(1,0);
subplot(333)
plot(Noise)
title("C")
ylim([-5,5]);
subplot(336)
hist(Noise)
subplot(339)
pwelch(Noise)
