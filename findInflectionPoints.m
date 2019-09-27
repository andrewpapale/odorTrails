% find inflection points in dnT-dT;

thr = 25;
nM = max(mouse);
nS = max(sess);
ninflection = nan(nR,nS);

for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS & ~edge & ~isstopped & V<100 & dT<25;
        
        if sum(k)>0
            
            if spotfound(iM,iS)==1
                k = k & t2s0;
            end
            
            dnT0 = dnT(k);
            dT0 = dT(k);
            
            
            putativeinflectionpts = find(diff(dnT0-dT0)==0);
            criteria = diff(putativeinflectionpts)<thr;
            while any(criteria)
                putativeinflectionpts(criteria)=[];
                criteria = diff(putativeinflectionpts)<thr;
            end
            ninflection(iM,iS) = length(putativeinflectionpts);
            
        end
    end
    disp(iM);
end
