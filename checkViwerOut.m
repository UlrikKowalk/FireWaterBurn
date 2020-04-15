clear;
close all;
clc;


sFileName = 'Viwer_Out.txt';

load(sFileName);
data = [Viwer_Out(1:2:end)', Viwer_Out(2:2:end)'];

plot(data);

audiowrite('ViwerOut.wav', data/max(max(abs(data))), 8000);