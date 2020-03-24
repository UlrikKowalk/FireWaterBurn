ccc;

sFilename = 'threshold.txt';
sFileAudio = 'media/0_recording_joined.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

stEntries = extractTextMarkers('media/0_recording_joined.txt');

nBlockSize = 512;
numTicks = 8;
nStep = 10;
nTheta = 360/nStep;

nSourcesMax = 10;
data = reshape(threshold, nSourcesMax, length(threshold)/nSourcesMax);

nBlocks = size(data, 2);
vTimeBlock = (1:nBlocks)*nBlockSize/nFs;
vTime = (1:length(vSignal))/nFs;
vDirections = linspace(0, 360-nTheta, nTheta);









hFig1 = figure();
hAx1 = subplot(311);
%imagesc(vTimeBlock, 0:nStep:(360-nTheta), data*nStep); 
%hAx1.YDir = 'reverse';

hold all;
for iSource = 1:2
   plot(vTimeBlock, (data(iSource, :)+1)*nStep, '.'); 
end
% plot(vTimeBlock, 5*mean(data+1)*nStep);
hold off;

hold on;
for iEntry = 1:length(stEntries)
    
    start = stEntries(iEntry).start;
    stop = stEntries(iEntry).end;
    source = stEntries(iEntry).source;
    
    switch(source)
        case 'U'
            h = fill([start, stop, stop, start], [0, 0, 360, 360], 'r');
        case 'B'
            h = fill([start, stop, stop, start], [0, 0, 360, 360], 'b');
        case 'UB'
            h = fill([start, stop, stop, start], [0, 0, 360, 360], 'm');
    end
    set(h,'facealpha',.5)
    
end
hold off;


xlim([vTimeBlock(1), vTimeBlock(end)]);
ylim([0, 360]);
xlabel('Time [s]');
ylabel('DOA [°]');
box on;


subplot(312);
plot(vTime, vSignal);
xlabel('Time [s]');
axis tight;


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
