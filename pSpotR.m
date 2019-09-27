% RingTransitionMatrix
% 2017-08-08 AndyP
% get transition matrix for found/not found trials

doTest = false;
nR = 100;
dR = 11.2; % cm/pixel

nM = max(mouse);
nS = max(sess);

iC=1;
iD=1;
iE=1;
iF=1;
foundR = [];
notR = [];
mfound = [];
mnot = [];
mrandf = [];
mrandn = [];
mnotf = [];
frandR = [];
nrandR = [];
ktime = 1500;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0 && t2s(iM,iS)~=-40
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = nx(k);
            y0 = ny(k);
            
            nearX0 = nearX(k);
            nearY0 = nearY(k);
            
            edge0 = edge;
            %disp(size(x0));
            
            kyes = spotfound(iM,iS)==1;% & conc0(iM,iS)==2;% & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            kno = spotfound(iM,iS)==0;% & conc0(iM,iS)==2;% & conc0(iM,iS)==iconc & bait0(iM,iS)==1;
            
            if kyes
                ktime = round(t2s(iM,iS)*50)+100;
            else % do not update ktime
            end
            
            %disp(min(length(x0),ktime));
            
            x0 = x0(1:min(length(x0),ktime));
            y0 = y0(1:min(length(y0),ktime));
            edge0 = edge0(1:min(length(x0),ktime));
            x0 = x0(~edge0);
            y0 = y0(~edge0);
            
            % rand data for shuffle control
            rng('shuffle');
            xR = 1+(1280-1).*rand(length(x0),1);
            yR = 1+(1024-1).*rand(length(y0),1);
            randdnT = sqrt((xR-xT).^2+(yR-yT).^2)./11.2;
            % create circles around spot
            if any(randdnT < 2)
                iE=iE+1;
            else
                iF=iF+1;
            end
            
            
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
                    foundR(iR,iC) = nansum(inR)./(area*length(x0));
                    mfound(1,iC)=iM;
                elseif kno
                    notR(iR,iD)=nansum(inR)./(area*length(x0));
                    mnot(1,iD)=iM;
                end
                
                if any(randdnT < 2)
                    frandR(iR,iE) = nansum(inRrand)./(area*length(xR));
                    mrandf(1,iE) = iM;
                else
                    nrandR(iR,iF) = nansum(inRrand)./(area*length(xR));
                    mrandn(1,iF) = iM;
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

pSR = nan(nR,nM);
pSrandR = nan(nR,nM);
hit = nan(nR,nM);
miss = nan(nR,nM);
hitR = nan(nM,1);
missR = nan(nM,1);
for iM=1:nM
    for iR=1:nR
        hit(iR,iM) = nansum(foundR(iR,mfound==iM))./sum(mfound==iM);
        miss(iR,iM) = nansum(notR(iR,mnot==iM))./sum(mnot==iM);
        hitR(iR,iM) = nansum(frandR(iR,mrandf==iM))./sum(mrandf==iM);
        missR(iR,iM) = nansum(nrandR(iR,mrandn==iM))./sum(mrandn==iM);
        pSR(iR,iM) = (hit(iR,iM)./(hit(iR,iM)+miss(iR,iM)));
        pSrandR(iR,iM) = hitR(iR,iM)./(hitR(iR,iM)+missR(iR,iM));
    end
end

F = figure(1); clf; hold on;
% plot(radii,pSR,'linewidth',2);
% plot(radii,pSrandR,'k','linewidth',2);
lineProps.col = {'r'};
mseb(radii,nanmean(pSR,2)',nanstd(pSR,[],2)',lineProps);
lineProps.col = {'k'};
mseb(radii,nanmean(pSrandR,2)',nanstd(pSrandR,[],2)',lineProps);
% legend('data','random');
xlabel('distance from spot (cm)','fontsize',21);
ylabel('probability of finding spot','fontsize',21);
title('Probability of finding Spot vs Distance from Spot','fontsize',24);
set(gca,'fontsize',18);



