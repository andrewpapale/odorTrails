% get distribution of spots
nM = max(mouse);
nS = max(sess);

xT = nan(nM,nS);
yT = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        kT = mouseT==iM & sessT==iS;
        
        xT(iM,iS) = nanmedian(xT1(kT));
        yT(iM,iS) = nanmedian(yT1(kT));
    end
end
        