% random position control

nShuf = 100;
nM = max(mouse);
nS = max(sess);

s2spotfound = nan(nM,nS,nShuf);
s2t2s = nan(nM,nS,nShuf);


for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            kT = mouseT==iM & sessT==iS;
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            if ~(xT<44.8 | xT>1280-44.8 | yT<44.8 | yT>1024-44.8) %#ok<OR2>
                k0 = k & ~edge & ~isstopped & V<100;
                for iShuf=1:nShuf
                    rng('shuffle');
                    xR = 44.8+(1280-44.8).*rand(length(k0),1);
                    yR = 44.8+(1024-44.8).*rand(length(k0),1);
                    randdnT = sqrt((xR-xT).^2+(yR-yT).^2)./11.2;
                    
                    ix=find(randdnT<2,1,'first');
                    if ~isempty(ix)
                        s2spotfound(iM,iS,iShuf)=1;
                        s2t2s(iM,iS,iShuf)=ix/50;
                    else
                        s2spotfound(iM,iS,iShuf)=0;
                    end
                end
            end
        end
        fprintf('Session %d/%d \n',iS,nS);
    end
    disp(iM);
end