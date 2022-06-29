
nM = 4;
nS = 25;

Type = nan(size(x));
Amp = nan(size(x));
Nstim = nan(size(x));
Nstim0 = [];
Amp0 = [];
Type0 = [];
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        nT = sT(iM,iS);
        if sum(k)>0
            kL = mouseL==iM & sessL==iS;
            fkL = find(kL);
            ampOut0 = ampOut(kL);
            timeOut0 = timOut(kL);
            t2s1 = t2s(mouseT==iM & sessT==iS);
            tff0 = tff(mouseT==iM & sessT==iS);
            for iT=1:nT
                lapt = find(t2s0==iT & k);
                Amp(lapt) = repmat(mode(ampOut0(timeOut0 > tff0(iT) & timeOut0 < t2s1(iT))),[length(lapt),1]);
                Type(lapt) = repmat(cond(iM,iS),[length(lapt),1]);
                Nstim0(iC)= sum(timeOut0 > tff0(iT) & timeOut0 < t2s1(iT));
                Amp0(iC) = mode(ampOut0(timeOut0 > tff0(iT) & timeOut0 < t2s1(iT)));
                Type0(iC) = cond(iM,iS);
                iC = iC+1;
            end
        end
    end
    disp(iM);
end