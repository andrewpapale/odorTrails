% get prevC

nM = 5;
nS = 105;
prevC = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        Tr0 = Tr(iM,iS);
        if Tr0~=1
            prevC(iM,iS)=C(iM,iS-1);
        end
    end
end
