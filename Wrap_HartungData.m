function [x,y,nx,ny,confn,confb,dx,dy,V,nV,dphi,dphib,sess,InfoMatrix] = Wrap_HartungData()
% 2019-04-17 AndyP
% Wrap_HartungData
% concatenate dataset from Jane Hartung

x = [];
y = [];
nx0 = [];
ny0 = [];
confn0 = [];
confb0 = [];
dx = [];
dy = [];
V = [];
nV = [];
dphi = [];
dphib = [];
sess = [];

InfoMatrix = {...
    '192'	'M'	'A'
    '193'	'M'	'A'
    '196'	'M'	'C'
    '197'	'M'	'C'
    '199'	'M'	'A'
    '202'	'M'	'C'
    '203'	'M'	'C'
    '207'	'M'	'C'
    '208'	'M'	'C'
    '213'	'M'	'A'
    '214'	'M'	'A'
    '194'	'F'	'C'
    '195'	'F'	'C'
    '200'	'F'	'C'
    '201'	'F'	'C'
    '204'	'F'	'A'
    '205'	'F'	'A'
    '206'	'F'	'A'
    '209'	'F'	'A'
    '210'	'F'	'A'
    '211'	'F'	'C'
    '212'	'F'	'C'
    };

D1 = dir('*_chunk1.mat');
D2 = dir('*_chunk2.mat');
Dmov = dir('*.avi');

nD = length(D1); % first chunk starts at frame 1
nC = length(D2);
nM = length(Dmov);

assert(nM~=0,'.avi movie files are missing from folder');
assert(nD==size(InfoMatrix,1),'data size mismatch');
assert(nM==size(InfoMatrix,1),'data size mismatch');


for iD=1:nD
    % read avi file to get frameRate
    nM1 = Dmov(iD).name;
    vr = VideoReader(nM1); %#ok<TNMLP>
    frameRate = vr.FrameRate;
    
    nD1 = D1(iD).name;
    nD1p = strsplit(nD1,'_');
    
    load(nD1,'bx','by','nx','ny','confb','confn');
    
    x = cat(1,x,bx);
    y = cat(1,y,by);
    nx0 = cat(1,nx0,nx); %#ok<*NODEF>
    ny0 = cat(1,ny0,ny);
    confb0 = cat(1,confb0,confb);
    confn0 = cat(1,confn0,confn);
    dx0 = foaw_diff(bx, 1/frameRate, frameRate, 0.2, 0.1);
    dy0 = foaw_diff(by, 1/frameRate, frameRate, 0.2, 0.1);
    dx = cat(1,dx,dx0');
    dy = cat(1,dy,dy0');
    V = cat(1,V,sqrt(dx0'.^2+dy0'.^2));
    dxn = foaw_diff(nx, 1/frameRate, frameRate, 0.2, 0.1);
    dyn = foaw_diff(ny, 1/frameRate, frameRate, 0.2, 0.1);
    nV = cat(1,nV,sqrt(dxn'.^2+dyn'.^2));
    dphi0 = zIdPhi1(dxn,dyn,1/frameRate,frameRate,0.2,0.1);
    dphib0 = zIdPhi1(dx0,dy0,1/frameRate,frameRate,0.2,0.1);
    dphi = cat(1,dphi,dphi0');
    dphib = cat(1,dphib,dphib0');
    sess = cat(1,sess,repmat(iD,[length(bx),1]));
    
    for iC=1:nC
        nC1 = D2(iC).name;
        nC1p = strsplit(nC1,'_');
        
        if strcmp(nD1p{1},nC1p{1}) && strcmp(nD1p{2},nC1p{1}) % add chunk 1 to chunk 2 for the same session
            
            load(nC1,'bx','by','nx','ny','confb','confn');
            
            x = cat(1,x,bx);
            y = cat(1,y,by);
            nx0 = cat(1,nx0,nx); %#ok<*NODEF>
            ny0 = cat(1,ny0,ny);
            confb0 = cat(1,confb0,confb);
            confn0 = cat(1,confn0,confn);
            dx0 = foaw_diff(bx, 1/frameRate, frameRate, 0.2, 0.1);
            dy0 = foaw_diff(by, 1/frameRate, frameRate, 0.2, 0.1);
            dx = cat(1,dx,dx0');
            dy = cat(1,dy,dy0');
            V = cat(1,V,sqrt(dx0'.^2+dy0'.^2));
            dxn = foaw_diff(nx, 1/frameRate, frameRate, 0.2, 0.1);
            dyn = foaw_diff(ny, 1/frameRate, frameRate, 0.2, 0.1);
            nV = cat(1,nV,sqrt(dxn'.^2+dyn'.^2));
            dphi0 = zIdPhi1(dxn,dyn,1/frameRate,frameRate,0.2,0.1);
            dphib0 = zIdPhi1(dx0,dy0,1/frameRate,frameRate,0.2,0.1);
            dphi = cat(1,dphi,dphi0');
            dphib = cat(1,dphib,dphib0');
            sess = cat(1,sess,repmat(iD,[length(bx),1]));
            
        end
    end
    disp(iD);
end

nx = nx0;
ny = ny0;
confb = confb0;
confn = confn0;

end
