nM = 4;
nS = 105;

for iM=2%1:nM
    for iS=12:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        
        x0 = x(k);
        y0 = y(k);
        nx0 = nx(k);
        ny0 = ny(k);
        
        F = figure(1); clf;
        plot(x0,y0,'.','markersize',1,'color',[0.5 0.5 0.5]);
        hold on;
        f2s0 = f2s(kT);
        f2f0 = f2f(kT);
        for iT=1:sum(kT)
            if iT==1
                r0 = (1:f2s0(iT));
            else
                r0 = (f2f0(iT-1):f2s0(iT));
            end
            P1 = plot(x0(r0),y0(r0),'k.-');
            P2 = plot(nx0(r0),ny0(r0),'b.','markersize',5);
            viscircles([cx0(iM,iS),cy0(iM,iS)],5.1744*1.5);
            title(sprintf('Wind? %d, Trial %d/%d, (Mouse %d, Sess %d)', Wind(iM,iS),iT,sum(kT),iM,iS),'fontsize',24);
            pause;
            delete(P1);
            delete(P2);
        end
        
    end
end