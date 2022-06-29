castingB = abs(Iorient)>45;
castingA = znidphi1>0.5;
for iM=1:5 
    for iS=1:20
        clf;
        kT=mouseT==iM & sessT==iS; 
        k=mouse==iM & sess==iS & ~edge; 
        if sum(k)>0
            subplot(1,2,1); hold on;
            
            scatter(nx(k),ny(k),10,znidphi1(k),'filled');
            plot(nx(k & castingA),ny(k & castingA),'rx');
            title('idphi');
            colorbar;
            caxis([-1 3]);
            axis off;
            plot(xT1(kT),yT1(kT),'g.','markersize',20); 
            
            subplot(1,2,2); hold on;
            
            scatter(nx(k),ny(k),10,abs(Iorient(k)),'filled');
            quiver(x(k),y(k),nx(k)-x(k),ny(k)-y(k));
            plot(nx(k & castingB),ny(k & castingB),'rx');
            title('Orient');
            colorbar;
            caxis([0 180]);
            axis off;
            plot(xT1(kT),yT1(kT),'g.','markersize',20); 
            
            pause;
        end
    end
end
