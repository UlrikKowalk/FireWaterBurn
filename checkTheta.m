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

numTicks = 8;
nStep = 10;
nTheta = 360/nStep;
a = reshape(threshold, nTheta, length(threshold)/nTheta);
nBlocks = size(a, 2);
vTimeBlock = (1:nBlocks)*nBlockSize/nFs;
vTime = (1:length(vSignal))/nFs;

stSources = struct('Theta', [],'Magnitude', [], 'State', []);

nPeaks = 5;
nDistanceLimit = 0.3;
vLoudness = zeros(nBlocks, 1);

hFig2 = figure();
hold on;

tic;

% nBlocks = 20;

for iBlock = 1 : nBlocks
   
    vBlock = abs(a(:, iBlock));
    vLoudness(iBlock) = rms(vBlock);
    vBlock = vBlock / max(vBlock);
    
    [mag, loc] = findpeaks(vBlock, 'NPeaks', nPeaks);
    
    if bQuad
        for iLoc = 1:length(loc)
            if (loc(iLoc) > 1)
                [p,y] = quadraticInterpolation(vBlock(loc(iLoc)-1), ...
                    vBlock(loc(iLoc)),vBlock(loc(iLoc)+1));
                mag(iLoc) = y;
                loc(iLoc) = loc(iLoc) + p;
            end
        end
    end
    
    if iBlock == 1
        
        for iLoc = 1:length(loc)
            stSources(iLoc).Theta = loc(iLoc);
            stSources(iLoc).Magnitude = mag(iLoc);
            stSources(iLoc).State = 2;
        end
        
    else
        
        stSourcesPrev = stSources;
        
        vErase = [];
        for iSource = 1 : length(stSources)
            if (stSources(iSource).State == 2)
                vErase(end+1) = iSource;
            end
        end
        stSources(vErase) = [];
        
       
        vLocationsPrev = zeros(size(stSourcesPrev, 2), 1);
        
        for iLoc = 1 : size(stSourcesPrev, 2)
            vLocationsPrev(iLoc) = stSourcesPrev(iLoc).Theta;
        end
    
        for iLoc = 1 : size(loc, 1)
           
            vDistance = abs(vLocationsPrev - loc(iLoc));
            [~, arg] = min(vDistance);
            if vDistance(arg) <= nDistanceLimit
                stSources(iLoc).Theta = loc(iLoc);
                stSources(iLoc).Magnitude = mag(iLoc);
                stSources(iLoc).State = 1;
            else 
                stSources(iLoc).Theta = loc(iLoc);
                stSources(iLoc).Magnitude = mag(iLoc);
                stSources(iLoc).State = 0;
            end
            
        end
        
        for iLocPrev = 1 : size(stSourcesPrev, 2)
           
            vDistance = abs(loc - stSourcesPrev(iLocPrev).Theta);
            [~, arg] = min(vDistance);
            if arg > nDistanceLimit
                stSources(end+1).Theta = stSourcesPrev(iLocPrev).Theta;
                stSources(end+1).Magnitude = stSourcesPrev(iLocPrev).Magnitude;
                stSources(end+1).State = 2;
            else 
                fprintf('EXCEPTION\n');
            end
            
        end
          
    end
    
    for iSource = 1 : size(stSources, 2)
        switch(stSources(iSource).State)
            case 0
                plot(iBlock, stSources(iSource).Theta, 'yx');
            case 1
                plot(iBlock, stSources(iSource).Theta, 'rx');
            case 2
                plot(iBlock, stSources(iSource).Theta, 'bx');
        end
    end
    
%     drawnow;
    
    
end
 
hold off;


% vDirections = linspace(0, 360-nTheta, nTheta);
% hFig1 = figure();
% hAx1 = subplot(311);
% imagesc(vTimeBlock, 0:nStep:(360-nTheta),  ((abs(a).^1))); 
% hAx1.YDir = 'reverse';
% 
% hold on;
% % for iBlock = 1:nBlocks
% %    plot(vTimeBlock(iBlock)*ones(nPeaks, 1), (mPeaks(iBlock, :))*nStep, 'r.'); 
% % end
% hold off;
% 
% 
% xlabel('Time [s]');
% ylabel('DOA [°]');



% subplot(312);
% plot(vTime, vSignal);
% xlabel('Time [s]');
% axis tight;


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