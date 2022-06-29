% get time on each arm

thr = 5;

onA = dTarm{1}<thr & dTarm{2}>thr & dTarm{3}>thr;
onB = dTarm{2}<thr & dTarm{1}>thr & dTarm{3}>thr;
onC = dTarm{3}<thr & dTarm{2}>thr & dTarm{1}>thr;

nM = max(mouse);
nS = max(sess);


tonA = nan(nM,nS);
tonB = nan(nM,nS);
tonC = nan(nM,nS);
lenA = nan(nM,nS);
lenB = nan(nM,nS);
lenC = nan(nM,nS);
thetaA = nan(nM,nS);
thetaB = nan(nM,nS);
mouse0 = nan(nM,nS);
sess0 = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k0)>0
            theta10 = theta1(k0);
            theta20 = theta2(k0);
            theta30 = theta3(k0);
            
            onA0 = onA(k0);
            onB0 = onB(k0);
            onC0 = onC(k0);
            
            tonA(iM,iS) = nansum(onA0)./50;
            tonB(iM,iS) = nansum(onB0)./50;
            tonC(iM,iS) = nansum(onC0)./50;
            
            lenA(iM,iS) = nansum(arm{1}(kT));
            lenB(iM,iS) = nansum(arm{2}(kT));
            lenC(iM,iS) = nansum(arm{3}(kT));
            
            thetaA(iM,iS) = mode(theta10);
            thetaB(iM,iS) = mode(theta20);
            thetaC(iM,iS) = mode(theta30);
            
            mouse0(iM,iS)= iM;
            sess0(iM,iS)=iS;
            
        end
    end
end


