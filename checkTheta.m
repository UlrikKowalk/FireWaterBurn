clear;
close all;
clc;

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
nDistanceLimit = 1;
vLoudness = zeros(nBlocks, 1);

hFig2 = figure();
hold on;

for iBlock = 1 : nBlocks
   
    vLoudness(iBlock) = rms(a(:, iBlock));
    a(:, iBlock) = a(:, iBlock) / max(abs(a(:, iBlock)));
    [mag, loc] = findpeaks(abs(a(:, iBlock)), 'NPeaks', nPeaks);
    
    if iBlock == 1
        
        for iLoc = 1:length(loc)
            stSources(iLoc).Theta = loc(iLoc);
            stSources(iLoc).Magnitude = mag(iLoc);
            stSources(iLoc).State = 2;
        end
        
    else
        
        stSourcesPrev = stSources;
       
        vLocationsPrev = zeros(size(stSourcesPrev, 2), 1);
        
        for iLoc = 1 : size(stSourcesPrev, 2)
            vLocationsPrev(iLoc) = stSourcesPrev(iLoc).Theta;
        end
    
        for iLoc = 1 : size(loc)
           
            vDistance = abs(vLocationsPrev - loc(iLoc))
            [~, arg] = min(vDistance);
            if arg <= nDistanceLimit
                stSources(iLoc).Theta = loc(iLoc);
                stSources(iLoc).Magnitude = mag(iLoc);
                stSources(iLoc).State = 1;
            else 
                stSources(iLoc).Theta = loc(iLoc);
                stSources(iLoc).Magnitude = mag(iLoc);
                stSources(iLoc).State = 0;
            end
            
        end
        
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
% ylabel('DOA [�]');



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


% colorbar;