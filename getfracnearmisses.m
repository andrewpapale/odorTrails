% fracnearmisses

hitThr = 2;
missThr = 10;
Nnearmiss = nan(nM,nS);
anynearmiss = nan(nM,nS);
nM = max(mouse);
nS = max(sess);
    for iM=1:nM
        for iS=1:nS
            k = mouse==iM & sess==iS & ~edge;
            
            if sum(k)>0
                
                t1 = find(dnT(k)>hitThr & dnT(k)<=missThr);
                thit = find(dnT(k)<=hitThr);
                
                if ~isempty(t1) && isempty(thit)
                    Nnearmiss(iM,iS)=length(t1);
                    anynearmiss(iM,iS)=1;
                else
                    Nnearmiss(iM,iS)=0;
                    anynearmiss(iM,iS)=0;
                end
            end
        end
        disp(iM);
    end