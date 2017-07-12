pos = dir('*_positions.mat');
stf = dir('*-StartFrame.mat');
trf = dir('*-Odor-Trail.mat');
nP = length(pos);
nS = length(stf);
nT = length(trf);

assert(nP==nS,'must be same number of position and start files');
assert(nP==nT,'must be same number of position and odor trail files');

homeDir = cd;

for iP=1:nP
    
    load(pos(iP).name,'position_results');
    load(stf(iP).name,'startFrame');
    load(trf(iP).name,'data');
    
    tempStr = strsplit(pos(iP).name,'_');
    movieName = strcat(tempStr{1},'_',tempStr{2},'.avi');
    
    cd ..
    
    movieFile = dir(movieName);
    disp(movieFile.name);
    v = VideoReader(movieFile.name); %#ok<TNMLP>
    video = read(v,[startFrame 1000+startFrame]); %#ok<VIDREAD>
    
    cd(homeDir);
    
    [x0,y0,nx0,ny0,V,nV] = Process_VT(position_results,startFrame);
    
    F = figure(2); clf;
    
    for iT=startFrame:1000+startFrame
        imagesc(squeeze(video(:,:,:,iT))); hold on;
        if iT>10
            P = plot(x0(iT-10:iT),y0(iT-10:iT),'r.','markersize',20);
            N = plot(nx0(iT-10:iT),ny0(iT-10:iT),'b.','markersize',20);
        else
            P = plot(x0(1:iT),y0(1:iT),'r.','markersize',20);
            N = plot(nx0(1:iT),ny0(1:iT),'b.','markersize',20);
        end
        
        title(mat2str(iT));
        drawnow;
        pause(0.01);
        delete(P);
        delete(N);
    end
    
    
    
    
end