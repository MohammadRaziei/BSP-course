disp('Question 1 is running ...');
%% global variables or settings
fs = 250;
t = ((1:2500)-1)/fs;
task1 = data{1}{4};
plot(t,task(1,:))
title(['C3',',',data{1}{1},',',data{1}{2},',',data{1}{3}]);xlabel('t(sec)');


% figure('units','normalized','outerposition',[0 0 1 1])
% annotation('line', 0.37*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);
% annotation('line', 0.65*[1,1], [0. 1], 'Color', 'r', 'LineWidth', 2);

%% A


