function TestOptoVideo(vr,dataFile)
% 2019-04-26 AndyP
% TestOptoVideo(fileName)

load(dataFile);

time1 = 0:1/vr.FrameRate:vr.Duration;
time1 = time1(1:end-1);

timeOut1 = interp1(time0,time1,timeOut0,'nearest');

saveStr0 = strsplit(vr.Name,'.');
saveStr = strcat(saveStr0{1},'-Annotated.avi');
newV = VideoWriter(saveStr);
open(newV);
nT = 9000;
for iT = 1:nT
    t0 = interp1(time1,frame0,vr.CurrentTime,'nearest');
    I0 = vr.readFrame();
    I0 = imadjust(I0);
    figure(1);
    imagesc(I0);
    colormap gray;
    hold on;
    viscircles([cx,cy],radiusAroundSpot*pixpercm);
    plot(x0(t0),y0(t0),'ro','markersize',5);
    plot(nx0(t0),ny0(t0),'gx','markersize',10);
    if t0 > 6
        plot(nx0(t0-5:t0-1),ny0(t0-5:t0-1),'b.','markersize',5);
    end
    if any(vr.CurrentTime-round(timeOut1,2)<=0.2 & vr.CurrentTime-round(timeOut1,2)>0)
        s1 = 'STIM';
    else
        s1 = '';
    end
    title(strcat(mat2str(round(time0(t0),2)),'...',s1),'fontsize',24);
    frame = getframe(gcf);
    writeVideo(newV,frame);
    clf;
end
close(newV);