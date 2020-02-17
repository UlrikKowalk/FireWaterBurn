clear;
close all;
clc;

sFilename = 'threshold.txt';
sFileAudio = 'media/45_8_Milbe.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);


nTheta = 18;
a = reshape(threshold, nTheta, length(threshold)/nTheta);

hFig1 = figure();
subplot(211);
imagesc(log10(a)); 
subplot(212);
plot(vSignal);
axis tight;
% colorbar;