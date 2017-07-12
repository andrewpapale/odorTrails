function Test_VT(x0,y0,nx0,ny0,arena_data)
% 2017-06-08 AndyP
% Test_VT(x0,y0,nx0,ny0,arena_data)
% to plot and test the processed body/nose coordinates

doVideo = false;

F = figure(2); clf;
clf; imagesc(arena_data.MedianImage); hold on;
map = colormap('gray');
if doVideo
    v = VideoWriter('test.avi');
    open(v);
end
for iT=1:length(x0);
    P = plot(x0(1:iT),y0(1:iT),'r.','markersize',20);
    N = plot(nx0(iT),ny0(iT),'b.','markersize',20);
    title(mat2str(iT));

    
    if doVideo
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    pause(0.01);
    delete(P);
    delete(N);
end

if doVideo
    close(v);
end
end

