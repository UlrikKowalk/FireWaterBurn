
nSensors = 15;


mCoordinates = [logspace(-2, -0.6, nSensors)', ...
            logspace(-0.6, -2, nSensors)', zeros(nSensors, 1)];
mC = mCoordinates - mCoordinates(8, :);

nAlpha = 45+90;
mR = [cosd(nAlpha), -sind(nAlpha); sind(nAlpha), cosd(nAlpha)];

mC = [(mR * mC(:, 1:2)')', zeros(nSensors, 1)];

for iCoord = 1:size(mC,1)
   
    fprintf(['this.coordinates[',num2str(iCoord-1),'] = new double[] {',...
        num2str(mC(iCoord,1)),'f,',num2str(mC(iCoord,2)),...
        'f,',num2str(mC(iCoord,3)),'f};\n']);
    
    
end