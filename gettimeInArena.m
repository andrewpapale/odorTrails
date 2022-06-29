% gettimeinarena

iC = 1;
time_in_Arena = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS & ~edge;
        if sum(k)>0 & t2s(iM,iS)~=-40
            time_in_Arena(iC,1)=nansum(k)/50;
            iC = iC+1;
        end
    end
end