clear;
close all;
clc;

bQuad = true;
bSave = false;

width = 16/2.54;     % Width in inches
height = 8/2.54;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 9;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize


sFilename = 'threshold_kal.txt';
sFileAudio = 'media/0_recording_joined.wav';

load(sFilename);

[vSignal, nFs] = audioread(sFileAudio);

nBlockSize = 512;
nSources = 5;

numTicks = 8;
nStep = 10;
nTheta = 360/nStep;

a = reshape(threshold_kal, nTheta+10+2*nSources, length(threshold_kal)/(nTheta+10+2*nSources));
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
   
end
 
hFig1 = figure();
hAx1 = subplot(311);
imagesc(vTimeBlock, (0:(nTheta-1))*nStep,  log10(a + (min(min(a)))*sign(min(min(a)))+1).^2); 

nMagMin = 1.2;

hold on;
for iBlock = 1:nBlocks
    
    vPks = mPeaks(1:nPeaks, iBlock);
    vPks(vPks == 0) = [];
    nPks = length(vPks);
    plot(vTimeBlock(iBlock)*ones(nPks, 1), vPks*nStep, 'rx'); 
   
%     for iPeak = 1:length(mMagsMat(1:nPeaks, iBlock))
%         if mMagsMat(iPeak, iBlock) > nMagMin
%             plot(vTimeBlock(iBlock)*ones(nPeaks, 1), (mPeaksMat(iPeak, iBlock)-1)*nStep, 'b.'); 
%         end
%     end
    
    for (iSource = 1:nSources)
        if (mKalman(iSource, iBlock) ~= 0)
            plot(vTimeBlock(iBlock), (mKalman(iSource, iBlock))*nStep, 'gd'); 
            %plot(vTimeBlock(iBlock), (mKalmanWeighted(iSource, iBlock))*nStep, 'y+');
            %plot(vTimeBlock(iBlock), (mKalOut(iBlock, iSource))*nStep, 'gd'); 
        end
    end
    
end
hold off;


% hAx1.YDir = 'reverse';
xlabel('Time [s]');
ylabel('DOA [°]');





subplot(312);
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


%% Plotting the restored image

sFileName = 'Viwer_Out.txt';

load(sFileName);
data = [Viwer_Out(1:2:end)', Viwer_Out(2:2:end)'];

subplot(313);
plot(data);
axis tight;


audiowrite('ViwerOut.wav', data/max(max(abs(data))), 8000);




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

set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties


if bSave
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*300, height*300]); %<- Set size
    set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

    % Here we preserve the size of the image when we save it.
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    
    %set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches');%,'PaperSize',[pos(3), pos(4)])

    print('Source_Tracking','-dpdf','-r300');
    close all;
end

fprintf('Done.\n');
toc;


% colorbar;