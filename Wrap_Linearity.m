%Wrap_Linearity


postSmoothing = 0.25; % s
%window = 0.5; % s
m = 50;
d = 1;
dtime = 1/30; % Hz
L = [];
nL = [];
positFiles = dir('*_positions.mat');
startFiles = dir('*-StartFrame.mat');

nD = length(positFiles);

for iD=1:nD
    
    fprintf('%s %d/%d \n',positFiles(iD).name,iD,nD);
    load(positFiles(iD).name,'position_results');
    load(startFiles(iD).name,'startFrame');
    
    [x0,y0,nx0,ny0] = Process_VT(position_results,startFrame);
    
    L0 = SplineCurvature(x0,y0);
    L = cat(1,L,L0');
    
    nL0 = SplineCurvature(nx0,ny0);
    nL = cat(1,nL,nL0');
    
end