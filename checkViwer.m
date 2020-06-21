clear;
close all;
clc;

sFilename = 'threshold.txt';
sFileAudio = 'media/Boomerang.wav';

% sFileAudio = 'media/0_recording_joined.wav';

mColor = getColor();

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

nMaxTheta = 180;
nMargin = 20;
nMaxSources = threshold(1);
nTheta = threshold(2);
threshold(1:2) = [];

data = reshape(threshold, nMargin+nMaxSources+1, length(threshold)/(nMargin+nMaxSources+1));

vNumSources = data(1, :);

vCandidates = data(2:nMargin+1, :);
vCandidates(vCandidates == 0) = NaN;
vKalmanPeaks = data(nMargin+2:end,:);
vKalmanPeaks(vKalmanPeaks == 0) = NaN;

vX = (1:length(vKalmanPeaks))*256/8000;

vC = vCandidates(1,:);
vC(isnan(vC)) = 0;
[p, v] = testkalPVA(vC, 0.032);


hFig1 = figure();
subplot(2,1,1);

yyaxis left;
hold on;
plot(vX, vCandidates' * (nMaxTheta / nTheta), '.', 'Color', mColor(1,:));
plot(vX, vKalmanPeaks' * (nMaxTheta / nTheta), '.', 'Color', mColor(5,:));
plot(vX, p * (nMaxTheta / nTheta), 'Color', mColor(2,:));
hold off;
ylabel('Source Directions');
xlabel('Time [s]');
ylim([0,nMaxTheta-1]);
xlim([vX(1), vX(end)]);
yyaxis right;
plot(vX, vNumSources, 'Color', mColor(2,:));
ylim([0, nMaxSources]);
set(gca, 'YTick',(0:nMaxSources));
ylabel('# Tracked Sources');
title('Tracked Sources');

subplot(2,1,2);

vTime = (1:length(vSignal))/nFs;
plot(vTime, vSignal);
xlabel('Time [s]');
axis tight;
title('Audio Signal');
ylabel('Signal');

if exist([sFileAudio(1:end-3),'txt'], 'File')
    stEntries = extractTextMarkers('media/0_recording_joined.txt');
    hold on;
    for iEntry = 1:length(stEntries)
        
        start = stEntries(iEntry).start;
        stop = stEntries(iEntry).end;
        source = stEntries(iEntry).source;
        
        switch(source)
            case 'U'
                h = fill([start, stop, stop, start], [-0.5, -0.5, 0.5, 0.5], 'r');
            case 'B'
                h = fill([start, stop, stop, start], [-0.5, -0.5, 0.5, 0.5], 'b');
            case 'UB'
                h = fill([start, stop, stop, start], [-0.5, -0.5, 0.5, 0.5], 'm');
        end
        set(h,'facealpha',.5)
        
    end
    hold off;
end




