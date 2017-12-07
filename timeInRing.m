% TimeInRing
% 2017-08-07 AndyP
% get the time in each ring around spots

nR = 31;
dR = 11.2; % cm/pixel
ktime = 1500; % frames

nM = max(mouse);
nS = max(sess);

R = nan(nR,nM,nS);
Rrat = nan(nR,nM,nS);
randR = nan(nR,nM,nS);
radii0 = nan(nR,nM,nS);
sess0 = nan(nR,nM,nS);
trial0 = nan(nR,nM,nS);
conc0 = nan(nR,nM,nS);
mouse0 = nan(nR,nM,nS);
ALorKP0 = nan(nR,nM,nS);
bait0 = nan(nR,nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = nx(k);
            y0 = ny(k);
            
            kyes = t2s(iM,iS)>0;% & conc0(iM,iS)==iconc & bait0(iM,iS)==0;
            kno = t2s(iM,iS)==-20;% & conc0(iM,iS)==iconc & bait0(iM,iS)==0;
            
            if kyes
                ktime = round(t2s(iM,iS)*50)+500;
            end
            
            
            x0 = x0(1:min(length(x0),ktime));
            y0 = y0(1:min(length(y0),ktime));
            
            % rand data for shuffle control
            xR = nanmin(x)+(nanmax(x)-nanmin(x)).*rand(length(x0),1);
            yR = nanmin(y)+(nanmax(y)-nanmin(y)).*rand(length(y0),1);
            
            % create circles around spot
            
            radius = dR.*(1:nR);
            
            for iR=1:nR
                dist = sqrt((x0-xT).^2+(y0-yT).^2);
                randdist = sqrt((xR-xT).^2+(yR-yT).^2);
                if iR>1
                    inR = dist <= radius(iR) & dist > radius(iR-1);
                    inRrand = randdist <= radius(iR) & randdist > radius(iR-1);
                    area = pi*(radius(iR).^2-radius(iR-1).^2);
                else
                    inR = dist <= radius(iR);
                    inRrand = randdist <= radius(iR);
                    area = pi*radius(iR).^2;
                end
                
                %area=1;
                R(iR,iM,iS) = nansum(inR)./(50*area);
                Rrat(iR,iM,iS) = nansum(inR)./(50*area*length(x0));
                
                randR(iR,iM,iS) = nansum(inRrand)./(50*area);
            end
            
            
            radii0(:,iM,iS)=radius./11.2;
            conc0(:,iM,iS)=repmat(mode(conc(k)),[nR,1]);
            sess0(:,iM,iS)=repmat(mode(sess(k)),[nR,1]);
            mouse0(:,iM,iS)=repmat(mode(mouse(k)),[nR,1]);
            ALorKP0(:,iM,iS)=repmat(mode(ALorKP(k)),[nR,1]);
            trial0(:,iM,iS)=repmat(mode(trial(k)),[nR,1]);
            bait0(:,iM,iS)=repmat(mode(bait(k)),[nR,1]);
        end
    end
end
radii=radius./11.2;