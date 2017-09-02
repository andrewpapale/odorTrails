function Test_VT(x0,y0,nx0,ny0,name,startFrame)
% 2017-06-08 AndyP
% Test_VT(x0,y0,nx0,ny0,arena_data)
% to plot and test the processed body/nose coordinates

doVideo = true;

F = figure(2); clf;
%clf; imagesc(arena_data.MedianImage); hold on;
if doVideo
    v1 = VideoWriter('test.avi');
    open(v1);
    
end

tempStr = strsplit(name,'_');
movieName = strcat(tempStr{1},'_',tempStr{2},'.avi');

movieFile = dir(movieName);
disp(movieFile.name);
v = VideoReader(movieFile.name);
video = read(v,[startFrame 1000+startFrame]); %#ok<VIDREAD>

for iT=startFrame:1000+startFrame;
    imagesc(squeeze(video(:,:,:,iT))); hold on;
    if iT>10
        P = plot(x0(iT-10:iT),y0(iT-10:iT),'r.','markersize',20);
        N = plot(nx0(iT-10:iT),ny0(iT-10:iT),'b.','markersize',20);
    else
        P = plot(x0(1:iT),y0(1:iT),'r.','markersize',20);
        N = plot(nx0(1:iT),ny0(1:iT),'b.','markersize',20);
    end
    title(mat2str(iT));
    
    
    if doVideo
        frame = getframe(gcf);
        writeVideo(v1,frame);
    end
    drawnow;
    pause(0.01);
    delete(P);
    delete(N);
end

if doVideo
    close(v1);
end
end

