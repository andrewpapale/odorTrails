% convert to polar <r,theta>
nM = max(mouse);
nS = max(sess);
r = nan(size(x));
theta = nan(size(x));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
           
            x0 = nx(k);
            y0 = ny(k);
            
            r(k) = sqrt((x0-xT(iM,iS)).^2+(y0-yT(iM,iS)).^2)./11.2;
            theta(k) = atan2((y0-yT(iM,iS)),(x0-xT(iM,iS)));
            
        end
    end
    disp(iM);
end