t2s0 = zeros(size(x));
nVtr = nan(size(nTr));
Type = nan(size(nTr));
condT = nan(size(nTr));
sPT = nan(size(nTr));
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = sessT == iS & mouseT == iM;
        kL = mouseL ==iM & sessL == iS;
        if sum(k)>0 && sum(kT)>0
            f2sN = f2s(kT);
            f2fN = f2f(kT);
            nTrN = nTr(kT);
            nV0 = nV(k);
            x0 = x(k);
            nT = sum(~isnan(nTrN));
            ampOut0 = ampOut(kL);
            fraOut0 = fraOut(kL);
            krem = edge(k) | nV(k) > 100 | bV(k) > 100;
            nV0(krem) = nan;
            tstart = find(k==1,1,'first');
            for iT=1:nT
                tend = f2sN(iT);
                if iT==1
                    tlocstart = 0;
                else
                    tlocstart = f2fN(iT-1);
                end
                if ~isnan(tend)
                    t2s0(tstart+tlocstart:tstart+tend)= repmat(nTrN(iT),[1,length(tstart+tlocstart:tstart+tend)]);
                    nVtr(iC) = nanmean(nV0(tlocstart+1:tend));
                    amp11 = ampOut0(fraOut0 > tlocstart & fraOut0 <= tend & ~isnan(x0(fraOut0)));
                    Type(iC) = mode(amp11);
                    condT(iC) = cond(iM,iS);
                    sPT(iC) = sP(iM,iS);
                    assert(tstart+tend <= find(k==1,1,'last')); 
                    iC = iC+1;
                end
            end
        end
    end
end