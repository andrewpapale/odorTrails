function znC = reZscore(nC,mouse,sess)

% re-zscore nC

nM = max(mouse);
nS = max(sess);

k0 = nC < 1E5;

znC = nan(size(nC));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            
            znC(k & k0) = RobustZ(nC(k & k0));
        end
    end
end
