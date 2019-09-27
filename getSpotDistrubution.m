% get distribution of spots
nM = max(mouse);
nS = max(sess);

xT = nan(nM,nS);
yT = nan(nM,nS);
oob = 0;
for iM=1:nM
    for iS=1:nS
        kT = mouseT==iM & sessT==iS;
        
        xT(iM,iS) = nanmedian(xT1(kT));
        yT(iM,iS) = nanmedian(yT1(kT));
        
        edge = xT(iM,iS) < 44.8 | xT(iM,iS) > (1280-44.8) | yT(iM,iS) < 44.8 | yT(iM,iS) > (1024-44.8);
        
        if edge
            oob = oob+1;
        end
    end
end
        