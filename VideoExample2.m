
V0 = VideoReader(V);
if ~isempty(V0)
    video = read(V0,[5800,6500]);
else
    error('unknown video');
end


nF = size(video,4);
dt = cat(1,0.02,diff(time0(5300:6000)));
dnT = sqrt((nx0-cx).^2+(ny0-cy).^2);

V1 = VideoWriter(strcat(V(1:19),'-annotated-1-Process_Hartung_Video'),'Uncompressed AVI');
set(V1,'FrameRate',round(1./nanmedian(diff(time0))));
open(V1);
for iF=1:nF
    F = figure(1);
    imagesc(squeeze(video(:,:,1,iF))); 
    colormap gray;
    axis xy;
    hold on;
    plot(x0(iF+5300),y0(iF+5300),'ro','markersize',10);
    plot(x0(iF+5300),y0(iF+5300),'rx','markersize',10);
    plot(nx0(iF+5300),ny0(iF+5300),'g.','markersize',10);
    %plot(tx0(iF+5300),ty0(iF+5300),'b.','markersize',10);
    viscircles([cx cy],pixpercm*radiusAroundSpot);
    text(500,50,mat2str(round((dnT(iF+5300)))));
    pause(dt(iF));
    frame = getframe(gcf);
    writeVideo(V1,frame);
    clf;
end
close(V1);

