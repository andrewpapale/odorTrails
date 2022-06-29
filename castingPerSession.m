% casts / session
sumCB = nan(nM,nS);
sumCA = nan(nM,nS);
meanCC = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS & ~edge & frame<t2s(iM,iS)*50+500;
        
        if sum(k)>0
            sumCA(iM,iS)=nansum(castingA(k));
            sumCB(iM,iS)=nansum(castingB(k));
            meanCC(iM,iS)=nanmean(zdphi1(k));
        end
    end
    disp(iM);
end
            