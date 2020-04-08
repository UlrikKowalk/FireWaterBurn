clear;
close all;
clc;

sFileSpec = 'Viwer_Spec.txt';
sFileSignal = 'Viwer_Signal.txt';
sFileSource = 'media/0_recording_joined.wav';

load(sFileSignal);
load(sFileSpec);

vSpec = Viwer_Spec(1:2:end) + 1i*Viwer_Spec(2:2:end);

vSpec_Full = [vSpec, conj(vSpec(end-1:-1:2))];
vSpec_IFFT = ifft(vSpec_Full);
% vSpec_IFFT = fftshift(vSpec_IFFT);

vSource = audioread(sFileSource);
vBlock = vSource(1:512, 1).*sqrt(hann(512));


hold all;
plot(vSpec_IFFT);
plot(Viwer_Signal);
plot(vBlock, '--');
hold off;