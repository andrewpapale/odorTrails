V = VideoReader('356074-2018-11-19-1.avi');
nF = length(bx);

%ix = interp1(time0,frame0,timeOut0,'nearest');

V1 = VideoWriter('newfile.avi','Grayscale AVI');
V1.FrameRate = 35;
open(V1);
iC=1;
for iF=1:nF
    F = figure(1); clf;
    vid = V.readFrame;
    vid = squeeze(vid(:,:,1));
    vid = imadjust(uint8(vid));
    imagesc(vid);
%     axis xy;
%     colormap gray;
%     saveas(F,'temp.jpg');
%     clf;
%     temp = imread('temp.jpg');
%     imagesc(temp);
    hold on;
    axis xy;
    colormap gray;
    axis equal;
    axis off;
    if bp(iF)>0.8
        plot(bx1(iF),by1(iF),'ro','markersize',10);
    end
    if np(iF)>0.8
        plot(nx1(iF),ny1(iF),'go','markersize',10);
    end
    if tbp(iF)>0.8
        plot(tbx(iF),tby(iF),'yo','markersize',10);
    end
    if confb(iF)>0.8
        plot(bx(iF),by(iF),'rx','markersize',10);
    end
    if confn(iF)>0.8
        plot(nx(iF),ny(iF),'gx','markersize',10);
    end
%     if blp(iF)>0.5
%         plot(blx(iF),bly(iF),'b.','markersize',20);
%     end
%     if brp(iF)>0.5
%         plot(brx(iF),bry(iF),'.','color',[0.1 0.25 0],'markersize',20);
%     end
%     if flp(iF)>0.5
%         plot(flx(iF),fly(iF),'.','color','c','markersize',20);
%     end
%     if frp(iF)>0.5
%         plot(frx(iF),fry(iF),'.','color','m','markersize',20);
%     end
    viscircles([cx,cy],radiusAroundSpot*pixpercm);
%     if any(ix==iF)
%         %scatter(x(iF),y(iF),amp0(iC)*50,255);
%         iC=iC+1;
%         %scatter(x0(iF),y0(iF),30,amp1(ix));
%         rectangle('Position',[50 400 50 50],'FaceColor','red');
%     end
    %caxis([2 4]);
    F1 = getframe(F);
    F1.cdata = F1.cdata(:,:,1);
    writeVideo(V1,F1);
end
close(V1);