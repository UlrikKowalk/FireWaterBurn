clear;
close all;
clc;

bQuad = true;

% STATES: 0: new, 1: alive, 2: dying

sFilename = 'threshold.txt';
sFileAudio = 'media/0_recording_joined.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

nBlockSize = 512;
nSources = 2;

numTicks = 8;
nStep = 10;
nTheta = 360/nStep;
a = reshape(threshold, nTheta+10+2*nSources, length(threshold)/(nTheta+10+2*nSources));
mPeaks = a(37:46,:);

nTau = nBlockSize/nFs;

mKalman = a(47:46+nSources, :);
mKalmanWeighted = a(47+nSources:46+2*nSources, :);

a(37:end,:) = [];

nBlocks = size(a, 2);
vTimeBlock = (1:nBlocks)*nBlockSize/nFs;
vTime = (1:length(vSignal))/nFs;

stSources = struct('Theta', [],'Magnitude', [], 'State', []);

nPeaks = 2;
nDistanceLimit = 0.3;
vLoudness = zeros(nBlocks, 1);

mPeaksMat = zeros(10, nBlocks);
mMagsMat = zeros(10, nBlocks);


aKal(1) = Kalman(nTau);
aKal(2) = Kalman(nTau);
mKalOut = zeros(nBlocks, nSources);

tic;

% nBlocks = 20;

for iBlock = 1 : nBlocks
   
    vLoudness(iBlock) = rms(a(:, iBlock));
    
%     a(:, iBlock) = a(:, iBlock) / max(a(:, iBlock));%vLoudness(iBlock).^1.4;
    a(:, iBlock) = a(:, iBlock) / vLoudness(iBlock);

    
    [mag, loc] = findpeaks(a(:, iBlock), 'NPeaks', nPeaks);
    
    if bQuad
        for iLoc = 1:length(loc)
            if (loc(iLoc) > 1)
                [p,y] = quadraticInterpolation(a(loc(iLoc)-1, iBlock), ...
                    a(loc(iLoc), iBlock),a(loc(iLoc)+1, iBlock));
                mMagsMat(iLoc, iBlock) = y;
                mPeaksMat(iLoc, iBlock) = loc(iLoc) + p;
            end
        end
    end
    
    
%     aKal(1).iterate(mPeaksMat(1));
%     aKal(2).iterate(mPeaksMat(2));
%     
%     mKalOut(iBlock, 1) = aKal(1).getData();
%     mKalOut(iBlock, 2) = aKal(2).getData();

%     
%     if iBlock == 1
%         
%         for iLoc = 1:length(loc)
%             stSources(iLoc).Theta = loc(iLoc);
%             stSources(iLoc).Magnitude = mag(iLoc);
%             stSources(iLoc).State = 2;
%         end
%         
%     else
%         
%         stSourcesPrev = stSources;
%         
%         vErase = [];
%         for iSource = 1 : length(stSources)
%             if (stSources(iSource).State == 2)
%                 vErase(end+1) = iSource;
%             end
%         end
%         stSources(vErase) = [];
%         
%        
%         vLocationsPrev = zeros(size(stSourcesPrev, 2), 1);
%         
%         for iLoc = 1 : size(stSourcesPrev, 2)
%             vLocationsPrev(iLoc) = stSourcesPrev(iLoc).Theta;
%         end
%     
%         for iLoc = 1 : size(loc, 1)
%            
%             vDistance = abs(vLocationsPrev - loc(iLoc));
%             [~, arg] = min(vDistance);
%             if vDistance(arg) <= nDistanceLimit
%                 stSources(iLoc).Theta = loc(iLoc);
%                 stSources(iLoc).Magnitude = mag(iLoc);
%                 stSources(iLoc).State = 1;
%             else 
%                 stSources(iLoc).Theta = loc(iLoc);
%                 stSources(iLoc).Magnitude = mag(iLoc);
%                 stSources(iLoc).State = 0;
%             end
%             
%         end
%         
%         for iLocPrev = 1 : size(stSourcesPrev, 2)
%            
%             vDistance = abs(loc - stSourcesPrev(iLocPrev).Theta);
%             [~, arg] = min(vDistance);
%             if arg > nDistanceLimit
%                 stSources(end+1).Theta = stSourcesPrev(iLocPrev).Theta;
%                 stSources(end+1).Magnitude = stSourcesPrev(iLocPrev).Magnitude;
%                 stSources(end+1).State = 2;
%             else 
%                 fprintf('EXCEPTION\n');
%             end
%             
%         end
%           
%     end
%     
%     for iSource = 1 : size(stSources, 2)
%         switch(stSources(iSource).State)
%             case 0
%                 plot(iBlock, stSources(iSource).Theta, 'yx');
%             case 1
%                 plot(iBlock, stSources(iSource).Theta, 'rx');
%             case 2
%                 plot(iBlock, stSources(iSource).Theta, 'bx');
%         end
%     end
    
%     drawnow;    
    
    
end
 



% vDirections = linspace(0, 360-nTheta, nTheta);
hFig1 = figure();
hAx1 = subplot(211);
% imagesc(vTimeBlock, 0:nStep:(360-nTheta),  (log10(((a)).^200))); 
imagesc(vTimeBlock, (0:(nTheta-1))*nStep,  log10(a + (min(min(a)))*sign(min(min(a)))+1).^2); 

nMagMin = 1.2;

hold on;
for iBlock = 1:nBlocks
    
    vPks = mPeaks(1:nPeaks, iBlock);
    vPks(vPks == 0) = [];
    nPks = length(vPks);
    plot(vTimeBlock(iBlock)*ones(nPks, 1), vPks*nStep, 'rx'); 
   
    for iPeak = 1:length(mMagsMat(1:nPeaks, iBlock))
        if mMagsMat(iPeak, iBlock) > nMagMin
            plot(vTimeBlock(iBlock)*ones(nPeaks, 1), (mPeaksMat(iPeak, iBlock)-1)*nStep, 'b.'); 
        end
    end
    
    for (iSource = 1:nSources)
        if (mKalman(iSource, iBlock) ~= 0)
            plot(vTimeBlock(iBlock), (mKalman(iSource, iBlock))*nStep, 'gd'); 
            plot(vTimeBlock(iBlock), (mKalmanWeighted(iSource, iBlock))*nStep, 'y+');
            %plot(vTimeBlock(iBlock), (mKalOut(iBlock, iSource))*nStep, 'gd'); 
        end
    end
    
end
hold off;


% hAx1.YDir = 'reverse';
xlabel('Time [s]');
ylabel('DOA [°]');





subplot(212);
plot(vTime, vSignal);
xlabel('Time [s]');
axis tight;

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


% 
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


fprintf('Done.\n');
toc;


% colorbar;