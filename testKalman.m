
clear;
close all;
clc;
 

nSpeed = 1;
nSamplerate = 100;
nLength = 2;
nStd = 0.5;
vTime = 0 : 1/nSamplerate : nLength;
vPositionReference = -10 + nSpeed*vTime;

nSamples = length(vTime);
vPositionMeasurement = vPositionReference + nSpeed*nStd^2 * randn(1, nSamples);
vVelocityMeasurement = [0, diff(vPositionMeasurement)] * nSamplerate;

oKalman = Kalman(1/nSamplerate);
vPositionEstimate = zeros(1, nSamples);
vVelocityEstimate = zeros(1, nSamples);

for iSample = 1 : nSamples
    
    nSample = vPositionMeasurement(iSample);
    oKalman.iterate(nSample)
    [vPositionEstimate(iSample), vVelocityEstimate(iSample)] = oKalman.getData();
    
end


hFig = figure();

hAx1 = subplot(2,1,1);

hold all;

plot(vTime, vPositionMeasurement, 'x');
plot(vTime, vPositionReference, '.');
plot(vTime, vPositionEstimate);

hold off;
box on;
hAx1.YLim = [vPositionReference(1) - 1, vPositionReference(end) + 1];

hAx2 = subplot(2,1,2);

hold all;

plot(vTime, vVelocityMeasurement, 'x');
plot(vTime, nSpeed * ones(1, nSamples), '.');
plot(vTime, vVelocityEstimate);

hold off;
box on;

hAx2.YLim = [nSpeed - 3*nStd*nSpeed, nSpeed + 3*nStd*nSpeed];