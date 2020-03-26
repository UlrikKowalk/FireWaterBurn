function nEstiPost = weightedPos(vCandidates, nTruePos, vInterval, nTheta)

% Cyclic probability 
nSpread = diff(vInterval);
vCandidates = [vCandidates-nSpread, vCandidates, vCandidates+nSpread];

vWeights = normpdf(vCandidates, nTruePos, nTheta);
vWeights = vWeights / sum(vWeights);
nEstiPost = sum(vCandidates.*vWeights);

end