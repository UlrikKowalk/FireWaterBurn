clear;
close all;
clc;

sFilename = 'threshold.txt';

load(sFilename);

nDirection = 45;

plot(threshold*360/2/pi);
hold on;
plot([1, length(threshold)], nDirection * [1, 1], 'k:');
hold off;

axis tight;
box on;

title('Evaluation of DU estimates');
ylabel('Angle [degrees]');