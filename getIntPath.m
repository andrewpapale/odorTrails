

intP = nan(size(t2s));
Idphi = nan(size(t2s));
iC = 1;
for iM=1:nM
    for iS=1:nS
        k0 = mouseT==iM & sessT==iS;
        if sum(k0)>0
            for iP=1:sum(k0)
                k = mouse==iM & sess==iS & ~edge & t2s0==iP;
                intP(iC) = nansum(abs(diff(nx(k)))+abs(diff(ny(k))))./5.174;
                Idphi(iC) = nansum(abs(ndphi(k))./nV(k));
            end
        end
        iC = iC+1;
    end
    disp(iM);
end