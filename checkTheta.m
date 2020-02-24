clear;
close all;
clc;

sFilename = 'threshold.txt';
sFileAudio = 'media/0_recording_joined.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

nBlockSize = 512;

numTicks = 8;
nStep = 10;
nTheta = 360/nStep;
a = reshape(threshold, nTheta, length(threshold)/nTheta);
nBlocks = size(a, 2);
vTimeBlock = (1:nBlocks)*nBlockSize/nFs;
vTime = (1:length(vSignal))/nFs;

cSources = {};

nPeaks = 5;
mPeaks = zeros(nBlocks, nPeaks);
mMags = zeros(nBlocks, nPeaks);
vLoudness = zeros(nBlocks, 1);

for iBlock = 1:nBlocks
   
    a(:, iBlock) = a(:, iBlock) / max(abs(a(:, iBlock)));
    vLocs = NaN * zeros(nPeaks, 1);
    vMags = NaN * zeros(nPeaks, 1);
    
    [mag, loc] = findpeaks(abs(a(:, iBlock)), 'NPeaks', nPeaks);
    
    vLocs(1:length(loc)) = loc;
    vMags(1:length(loc)) = mag;
    
    mPeaks(iBlock, :) = vLocs;
    mMags(iBlock, :) = vMags;
    vLoudness(iBlock) = rms(a(:, iBlock));
    
 end


vDirections = linspace(0, 360-nTheta, nTheta);
hFig1 = figure();
hAx1 = subplot(311);
imagesc(vTimeBlock, 0:nStep:(360-nTheta),  ((abs(a).^1))); 
hAx1.YDir = 'reverse';

hold on;
for iBlock = 1:nBlocks
   plot(vTimeBlock(iBlock)*ones(nPeaks, 1), (mPeaks(iBlock, :))*nStep, 'r.'); 
end
hold off;


xlabel('Time [s]');
ylabel('DOA [°]');



subplot(312);
plot(vTime, vSignal);
xlabel('Time [s]');
axis tight;


% mMags(isnan(mMags)) = 0;
% 
% subplot(313);
% hold all; 
% for iPeak = 1:nPeaks
%     plot(vTimeBlock, mMags(:, iPeak) / sum(mMags, 2));
% end
% box on;
% xlabel('Time [s]');
% ylabel('Magnitude (norm)');
% axis tight;
% 
% hold off;




% colorbar;