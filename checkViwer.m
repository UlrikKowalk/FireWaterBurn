clear;
close all;
clc;

sFilename = 'threshold.txt';
sFileAudio = 'media/0_recording_joined.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

nSources = threshold(1);
nTheta = threshold(2);
threshold(1:2) = [];

data = reshape(threshold, length(threshold)/nSources, nSources);
hold on
plot(data(data > 0) * (360 / nTheta), '.');