% RingTransitionMatrix
% 2017-08-08 AndyP
% get transition matrix for found/not found trials

doTest = false;
nR = 61;
dR = 11.2*2; % cm/pixel

nM = max(mouse);
nS = max(sess);

iC=1;
iD=1;
iE=1;
iF = 1;
foundR = [];
notR = [];
ktime = round(nanmedian(t2s(t2s>0))*50);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = nx(k);
            y0 = ny(k);
            
            %disp(size(x0));
            
            kyes = t2s(iM,iS)>0;% & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            kno = t2s(iM,iS)==-20;% & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            
            if kyes
                ktime = round(t2s(iM,iS)*50)+500;
            else % do not update ktime
            end
            
            %disp(min(length(x0),ktime));
            
            x0 = x0(1:min(length(x0),ktime));
            y0 = y0(1:min(length(y0),ktime));
            
            % rand data for shuffle control
            rng('shuffle');
            xR = 1+(1280-1).*rand(length(x0),1);
            yR = 1+(1024-1).*rand(length(y0),1);
            
            % create circles around spot
            
            radius = linspace(dR,dR*nR,nR);
            
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
                    foundR(iR,iC) = nansum(inR);
                    
                elseif kno
                    notR(iR,iD)=nansum(inR);
                    
                end
                
                
                if rand(1)>0.5
                    frandR(iR,iE) = nansum(inRrand);
                    iE=iE+1;
                else
                    nrandR(iR,iF) = nansum(inRrand);
                    iF = iF+1;
                end
                
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
        end
        
    end
    disp(iM);
end
radii=radius./(11.2);

pSR = nan(nR,1);
pSrandR = nan(nR,1);
for iR=1:nR
    hit = nansum(foundR(iR,:));
    miss = nansum(notR(iR,:));
    hitR = nansum(frandR(iR,:));
    missR = nansum(nrandR(iR,:));
    pSR(iR,1) = (hit./(hit+miss));
    pSrandR(iR,1) = hitR./(hitR+missR);
end

F = figure(1); clf; hold on;
plot(radii,pSR,'linewidth',2);
plot(radii,pSrandR,'k','linewidth',2);
legend('data','random');
xlabel('distance from spot (cm)','fontsize',21);
ylabel('probability of finding spot','fontsize',21);
title('Probability of finding Spot vs Distance from Spot','fontsize',24);
set(gca,'fontsize',18);



