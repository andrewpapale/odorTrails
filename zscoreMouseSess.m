function zOut = zscoreMouseSess(zIn,mouse,sess)

nM = max(mouse);
nS = max(sess);


ztemp = nan(size(zIn));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            ztemp(k)=nanzscore(zIn(k));
        end
    end
end

zOut = nan(size(zIn));
for iM=1:nM
    k = mouse==iM;
    if sum(k)>0
        zOut(k) = nanzscore(ztemp(k));
    end
end



end