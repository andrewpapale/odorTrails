% plotExampleOptoNoise
% 2019-07-30 AndyP



nM = 3;
nS = 40;
iM = randi(5);
iS = randi(40);

for ix=1:20
        iL = randi(sT(iM,iS));
        kL = t2s0==iL;
        k = mouse==iM & sess==iS & t2s0==iL;
        kS = mouseT==iM & sessT==iS;
        clf;
        if sum(k)>0 % then plot
            tff0 = tff(kS);
            tff0 = tff0;
            t2s2 = t2s(kS);
            t2s2 = t2s(kS);
            kL = mouseL==iM & sessL==iS;
            timeOut0 = timOut(kL);
            timeOut0 = timeOut0(timeOut0 > tff0(iL) & timeOut0 < t2s2(iL));
            x0 = x(k);
            y0 = y(k);
            bV0 = bV(k);
            plot(x0,y0,'k.');
            hold on;
            scatter(x0,y0,20,bV0,'filled');
            
            viscircles([cx0(iM,iS),cy0(iM,iS)],2*5.174);
        end
end