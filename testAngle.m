v = VideoReader('163241_2017-03-21-102037-0000.avi');
video = read(v,[279 3000]);

kT = mouseT==1 & sessT==1; k = mouse==1 & sess==1;
x0 = x(k);
y0 = y(k);
I0 = Iorient(k);
nx0 = nx(k);
ny0 = ny(k);

iFr = 1;
while iFr > 0
    iFr = input('Enter frame: ');
    temp = gpuArray(squeeze(video(:,:,:,iFr)));
    temp = imadjust(temp,[0; 0.25],[]);
    imagesc(temp);
    hold on;
    axis xy;
    plot(xT1(kT),yT1(kT),'r.','markersize',20);
    plot(x0(iFr),y0(iFr),'r.');
    plot(nx0(iFr),ny0(iFr),'g.');
    xI = cos(I0(iFr)*pi/180);
    yI = sin(I0(iFr)*pi/180);
    quiver(x0(iFr),y0(iFr),xI*100,yI*100,'color','b','linewidth',2);
end