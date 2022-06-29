% get prior trial as baited=1 or unbaited=0

nM = size(trial0,1);
nS = size(trial0,2);
nT = length(x);
prevT = nan(nT,1);
prevT0 = nan(nM,nS);
for iM=1:nM
    for iS=2:nS
        k = mouse==iM & sess==iS;
        if trial0(iM,iS)~=1;
            if bait0(iM,iS-1)==0
                prevT0(iM,iS)=0;
            elseif bait0(iM,iS-1)==1
                prevT0(iM,iS)=1;
            end
            prevT(k==1)=repmat(prevT0(iM,iS),[sum(k),1]);
        end
    end
end