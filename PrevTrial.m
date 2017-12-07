% get PrevTrial bait/unbait and concentration



nM = max(mouse);
nS = max(sess);
prevbait = nan(nM,nS);
prevconc = nan(nM,nS);
prevtrial = nan(nM,nS);
for iM=1:nM
    for iS=2:nS
        k = mouse==iM & sess==iS;
        kprev = mouse==iM & sess==iS-1;
        
        if sum(k)>0 & mode(trial(k==1))~=1
            prevbait(iM,iS)=mode(bait(kprev==1));
            prevconc(iM,iS)=mode(conc(kprev==1));
            prevtrial(iM,iS)=mode(trial(kprev==1));
        end
    end
    disp(iM);
end

            