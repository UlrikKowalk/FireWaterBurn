clear; 
close all;
clc;


sFolder = [pwd, filesep, 'media'];
[cFiles, sPath] = uigetfile('*.wav', 'Multiselect', 'On');
nFiles = length(cFiles);

sFileOut = [sFolder, filesep, cFiles{1}(1:end-4), '_joined.wav'];

if (nFiles == 1 || nFiles == 0)
    return;
end

for iFile = 1:length(cFiles)
    
    [mTmp, nFs] = audioread(fullfile(sPath, cFiles{iFile}));
    
    if (iFile == 1)
        mSignal = mTmp;
    else
        mSignal(end+1:end+size(mTmp, 1), :) = mTmp;
    end
    
end

audiowrite(sFileOut, mSignal, nFs);

fprintf("%d files successfully joined to: %s.\n", nFiles, sFileOut);