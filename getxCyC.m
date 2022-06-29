% get xC,yC for each session

positFiles = dir('*_positions.mat');
compFiles = dir('*-ConnectingPoint.mat');
startFiles = dir('*-StartFrame.mat');

nD = length(positFiles);

xC1 = [];
yC1 = [];
for iD=1:nD
    fprintf('%s %d/%d \n',positFiles(iD).name,iD,nD);
    load(positFiles(iD).name,'position_results');
    load(compFiles(iD).name,'xC','yC');
    load(startFiles(iD).name,'startFrame');
    
    [x0,y0,nx0,ny0] = Process_VT(position_results,startFrame);
    
    nP = length(x0);
    
    xC1 = cat(1,xC1,repmat(xC,[length(x0),1]));
    yC1 = cat(1,yC1,repmat(yC,[length(x0),1]));
end

