%shuffle spots

shufdnT = [];

nM = max(mouse);
nS = max(sess);

xT = [];
yT = [];

for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            kT = mouseT==iM & sessT==iS;
            
            % get all spots
            xT = cat(1,xT,nanmedian(xT1(kT)));
            yT = cat(1,yT,nanmedian(yT1(kT)));

        end
    end
    disp(iM);
end

% shuffle spots
rng('shuffle');
idx = randperm(length(xT));
xT2 = xT(idx);
yT2 = yT(idx);
iC = 1;
for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            
            nx0 = nx(k);
            ny0 = ny(k);
            
            xT0 = xT(iC);
            yT0 = yT(iC);
            iC = iC+1;
            
            % calculate distance from trail
            dnT0 = sqrt((nx0-xT0).^2+(ny0-yT0).^2);
            
            shufdnT = cat(1,shufdnT,dnT0);
        end
    end
    disp(iM);
end
            
            