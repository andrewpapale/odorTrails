% TimeInRing
% 2017-08-07 AndyP 
% get the time in each ring around spots

nR = 21;
dR = 11.2; % cm/pixel

nM = max(mouse);
nS = max(sess);

R = nan(nR,nM,nS);
randR = nan(nR,nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        xT = nanmedian(xT1(kT));
        yT = nanmedian(yT1(kT));
        
        x0 = x(k);
        y0 = y(k);
        
        % rand data for shuffle control
        xR = nanmin(x)+(nanmax(x)-nanmin(x)).*rand(length(x0),1);
        yR = nanmin(y)+(nanmax(y)-nanmin(y)).*rand(length(y0),1);
        
        % create circles around spot
        
        radius = dR.*(1:nR);
        
        for iR=1:nR-1
            dist = sqrt((x0-xT).^2+(y0-yT).^2);
            inR = dist > radius(iR) & dist < radius(iR+1);
            %area = pi*(radius(iR+1).^2-radius(iR).^2);
            area=1;
            R(iR,iM,iS) = nansum(inR)./(50*area);
            randdist = sqrt((xR-xT).^2+(yR-yT).^2);
            inRrand = randdist > radius(iR) & randdist < radius(iR+1);
            randR(iR,iM,iS) = nansum(inRrand)./(50*area);
        end
    end
end
radii=radius./11.2;