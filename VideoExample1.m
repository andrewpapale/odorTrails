
V0 = VideoReader(V);
clear video;
if ~isempty(V0)
    video = read(V0,[5800,8500]);
else
    error('unknown video');
end

k = 5800:8500;
nF = size(video,4);
%dt = cat(1,0.02,diff(time0(500:1500)));
dt = 1/V0.FrameRate;
%dnT0 = dnT(k);
x1 = x0(k);
y1 = y0(k);
nx1 = tnx(k);
ny1 = tny(k);
%cx0 = mode(cx(k));
%cy0 = mode(cy(k));

%V1 = VideoWriter(strcat(V(1:19),'-annotated-1-optimouse'),'Uncompressed AVI');
V1 = VideoWriter(strcat(V(1:11),'-annotated-1-Process_Hartung_Video'),'Uncompressed AVI');
%set(V1,'FrameRate',round(1./nanmedian(diff(time0))));
set(V1,'FrameRate',V0.FrameRate);
open(V1);
for iF=1:nF
    F = figure(1);
    imagesc(imadjust(squeeze(video(:,:,1,iF)))); 
    colormap gray;
    axis xy;
    hold on;
    plot(x1(iF),y1(iF),'ro','markersize',10);
    plot(x1(iF),y1(iF),'rx','markersize',10);
    plot(nx1(iF),ny1(iF),'g.','markersize',10);
    %plot(tx0(iF+500),ty0(iF+500),'b.','markersize',10);
    %viscircles([cx0 cy0],11.2*2);
    %text(500,50,mat2str(round((dnT0(iF)))));
    %pause(dt(iF));
    pause(dt);
    frame = getframe(gcf);
    writeVideo(V1,frame);
    clf;
end
close(V1);

