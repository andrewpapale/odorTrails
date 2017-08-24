% RingTransitionMatrix
% 2017-08-08 AndyP
% get transition matrix for found/not found trials

doTest = false;
nR = 31;
dR = 11.2; % cm/pixel
ktime = 1500; % frames

nM = max(mouse);
nS = max(sess);

iC=1;
iD=1;
iE=1;
foundR = [];
notR = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0;
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = nx(k);
            y0 = ny(k);
            
            %disp(size(x0));
            
            kyes = t2s(iM,iS)>0 & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            kno = t2s(iM,iS)==-20 & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            
            if kyes
                ktime = round(t2s(iM,iS)*50)+500;
            end
            
            %disp(min(length(x0),ktime));
            
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
                
                
                if kyes
                    foundR(iR,iC) = nansum(inR)./(50*area);
                elseif kno
                    notR(iR,iD)=nansum(inR)./(50*area);
                end
                randR(iR,iE) = nansum(inRrand)./(50*area);
                
                if doTest
                    F = figure(1); clf;
                    plot(x0,y0,'k.');
                    hold on;
                    plot(xT,yT,'r.','markersize',20);
                    circle([xT,yT],radius(iR),1000);
                    if iR>1
                        circle([xT,yT],radius(iR-1),1000);
                    end
                    plot(x0(inR),y0(inR),'r.','markersize',5);
                    pause;
                end
                
                
                
            end
            if kyes
                iC=iC+1;
            elseif kno
                iD=iD+1;
            end
            iE=iE+1;
        end
        
    end
end
radii=radius./11.2;

for iR=1:nR
    for iQ=1:nR
        Cfound(iR,iQ)=corr2(foundR(iR,:),foundR(iQ,:));
        %Cnot(iR,iQ) = corr2(notR(iR,:),notR(iQ,:));
        %Crand(iR,iQ)= corr2(randR(iR,:),randR(iQ,:));
    end
end
