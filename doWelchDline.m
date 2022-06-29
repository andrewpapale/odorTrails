% doWelchDline

nM = max(mouse);
nS = max(sess);

for iM=1:nM;
    for iS=1:nS;
        k = mouse==iM & sess==iS;
        
        if sum(k)>0 % GO
            
            k = k & ~edge
            Dline0 = Dline(k);