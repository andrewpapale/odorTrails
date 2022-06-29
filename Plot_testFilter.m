% test filter

nT = length(nx0);
dnew = cat(1,nan,sqrt(diff(tnx).^2+diff(tny).^2));
figure(1); clf;
plot(x0,y0,'.','markersize',1,'color',[0.5 0.5 0.5]);
hold on;
for iT=1:nT
    if iT>5
        P1 = plot(x0(iT-5:iT),y0(iT-5:iT),'ko-','markersize',10);
        P5 = plot(nx0(iT-1:iT),ny0(iT-1:iT),'gx-','markersize',5);
        P6 = plot(tnx(iT-1:iT),tny(iT-1:iT),'bs-','markersize',5);
    else
        P1 = [];
        P5 = [];
        P6 = [];
    end
    P2 = plot(x0(iT),y0(iT),'ro','markersize',10);
    P3 = plot(nx0(iT),ny0(iT),'gx','markersize',5,'linewidth',2);
    P4 = plot(tnx(iT),tny(iT),'bs','markersize',5,'linewidth',2);
    
    title(sprintf('Frame %d, x(t)=%0.2f (%0.2f)',frame0(iT),x0(iT),time0(iT)),'fontsize',24);
    
    
    if dnew(iT) < 10 || isnan(dnew(iT))
    else
        pause;
    end
    delete(P1);
    delete(P2);
    delete(P3);
    delete(P4);
    delete(P5);
    delete(P6);
end