% 2018-09-20 AndyP
% script to check t2s

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if t2s(iM,iS)>0 && t2s(iM,iS)<3
            x0 = x(k); y0 = y(k); 
            F = figure(1); clf; 
            scatter(x(k),y(k),20,frame(k)); hold on; 
            plot(nx(k),ny(k),'k.-'); 
            viscircles([cx0(iM,iS),cy0(iM,iS)],11.2*2); 
            X = find(dnT(k)<2,1,'first'); 
            disp(X/50); 
            plot(x0(X),y0(X),'ro','markersize',10); 
            plot(initX(iM,iS),initY(iM,iS),'bx','markersize',25); 
            pause; 
        end
    end
end