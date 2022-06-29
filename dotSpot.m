% get heading relative to trail

nM = max(mouse);
nS = max(sess);

tV = nan(size(dnT));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)~=0
            kT = mouseT==iM & sessT==iS;
            mxT = nanmedian(xT1(kT));
            myT = nanmedian(yT1(kT));
            
            nx0 = nx(k);
            ny0 = ny(k);
            x0 = x(k);
            y0 = y(k);
            
            dx = (nx0-x0);
            dy = (ny0-y0);
            mag = sqrt(dx.^2+dy.^2);
            vec1 = [dx./mag,dy./mag];
            
            dTx = (mxT-x0);
            dTy = (myT-y0);
            magT = sqrt(dTx.^2+dTy.^2);
            vec2 = [dTx./magT,dTy./magT];
            
            tV(k) = arccos(dot(vec1,vec2,2));
        end
    end
end