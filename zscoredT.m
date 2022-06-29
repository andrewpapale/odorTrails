nM = max(mouse);
nS = max(sess);

zdT = nan(size(dT));
for iM=1:nM
        k = mouse==iM;
        zdT(k)=nanzscore(dT(k));     
end