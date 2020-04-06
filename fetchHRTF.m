
clear;
close all;
clc;

tic;

nElevation = 0;
nDistance = 1.25;
nFFTLen = 512;

sFileName = 'CH07IK25.mat';
sFileName_Out = ['CH07IK25_', num2str(nElevation)];
load(sFileName);

nResAzimuth = 5;
nResElevation = 15;
nSets = 360/nResAzimuth;

mHRTF_L_real = zeros(nSets, nFFTLen/2+1);
mHRTF_L_imag = zeros(nSets, nFFTLen/2+1);

mHRTF_R_real = zeros(nSets, nFFTLen/2+1);
mHRTF_R_imag = zeros(nSets, nFFTLen/2+1);

mPositions = zeros(nSets, 3);
iIdx = 1;

for iPos = 1:length(SourcePositions)
   
    if (SourcePositions(iPos, [2,3]) == [nElevation, nDistance])
        tmp = fft(hrir_l(iPos, :));
        mHRTF_L(iIdx, :) = tmp(1:nFFTLen/2+1);
        mHRTF_L_real(iIdx, :) = real(tmp(1:nFFTLen/2+1));
        mHRTF_L_imag(iIdx, :) = imag(tmp(1:nFFTLen/2+1));
        
        tmp = fft(hrir_r(iPos, :));
        mHRTF_R(iIdx, :) = tmp(1:nFFTLen/2+1);
        mHRTF_R_real(iIdx, :) = real(tmp(1:nFFTLen/2+1));
        mHRTF_R_imag(iIdx, :) = imag(tmp(1:nFFTLen/2+1));
        
        mPositions(iIdx, :) = SourcePositions(iPos, :);
        iIdx = iIdx + 1;
    end
   
end

mPositions = flipud(circshift(mPositions, [-19, 0]));
mPositions(:, 2:3) = [];

mHRTF_L_real = circshift(mHRTF_L_real, [-19, 0, 0]);
mHRTF_L_imag = circshift(mHRTF_L_imag, [-19, 0, 0]);
mHRTF_R_real = circshift(mHRTF_R_real, [-19, 0, 0]);
mHRTF_R_imag = circshift(mHRTF_R_imag, [-19, 0, 0]);

vDirections = mPositions;

save(['hrir', filesep, sFileName_Out, '.mat'], 'mHRTF_L', 'mHRTF_R');
save(['hrir', filesep, 'directions.mat'], 'vDirections');
% save(['hrir', filesep, sFileName_Out, '.mat'], 'mHRTF_L_real', 'mHRTF_L_imag', 'mHRTF_R_real', 'mHRTF_R_imag');%, 'mPositions');
% save(['hrir', filesep, sFileName_Out, 'L.mat'], 'mHRTF_L');%, 'mPositions');
% save(['hrir', filesep, sFileName_Out, 'R.mat'], 'mHRTF_R');%, 'mPositions');

toc;