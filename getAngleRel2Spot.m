% get AngleRel2Spot

nM = max(mouse);
nS = max(sess);

nP = length(x);

angleRel2Spot = nan(nP,1);
for iM=1:nM
    disp(iM);
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            x1 = nx(k)-x(k);
            y1 = ny(k)-y(k);
            mag1 = sqrt(x1.^2+y1.^2);
            vec1 = [x1./mag1 y1./mag1];
            
            xT2 = nanmedian(xT1(kT));
            yT2 = nanmedian(yT1(kT));
            
            x2 = xT2-x(k);
            y2 = yT2-y(k);
            mag2 = sqrt(x2.^2+y2.^2);
            vec2 = [x2./mag2 y2./mag2];
            
            angleRel2Spot(k) = acos(dot(vec1,vec2,2))*180/pi;
            
        end
    end
end
