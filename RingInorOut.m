% RingInorOut
% 2017-08-10 AndyP
% get probability of going to nearer or father ring

doTest = false;
nR = 31;
dR = 11.2; % cm/pixel
ktime = 1500; % frames

nM = max(mouse);
nS = max(sess);


pIn = nan(nR,nM,nS);
pOut = nan(nR,nM,nS);
pInRand = nan(nR,nM,nS);
pOutRand = nan(nR,nM,nS);
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
            
            kyes = t2s(iM,iS)>0;% & conc0(iM,iS)==iconc & bait0(iM,iS)==0;
            kno = t2s(iM,iS)==-20;% & conc0(iM,iS)==iconc & bait0(iM,iS)==0;
            
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
            
            indx = repmat((1:length(x0))',[1,nR]);
            inR = nan(length(x0),1);
            inRrand = nan(length(x0),1);
            for iR=1:nR
                dist = sqrt((x0-xT).^2+(y0-yT).^2);
                randdist = sqrt((xR-xT).^2+(yR-yT).^2);
                
                if iR>1
                    inR(:,iR) = dist <= radius(iR) & dist > radius(iR-1);
                    inRrand(:,iR) = randdist <= radius(iR) & randdist > radius(iR-1);
                    %area = pi*(radius(iR).^2-radius(iR-1).^2);
                else
                    inR(:,iR) = dist <= radius(iR);
                    inRrand(:,iR) = randdist <= radius(iR);
                    %area = pi*radius(iR).^2;
                end
                
                %area=1;
                
                
                %                 if kyes
                %                     foundR(iR,iC) = nansum(inR)./(50*area);
                %                 elseif kno
                %                     notR(iR,iD)=nansum(inR)./(50*area);
                %                 end
                %                 randR(iR,iE) = nansum(inRrand)./(50*area);
                
                if doTest
                    F = figure(1); clf;
                    plot(x0,y0,'k.');
                    hold on;
                    plot(xT,yT,'r.','markersize',20);
                    circle([xT,yT],radius(iR),1000);
                    if iR>1
                        circle([xT,yT],radius(iR-1),1000);
                    end
                    plot(x0(inR{iR}),y0(inR{iR}),'r.','markersize',5);
                    pause;
                end
            end
            
            %keyboard;
            [indx0,ring]=find(inR==1);
            [~,indx1] = sort(indx0);
            ring = ring(indx1);
            
            diffring = cat(1,nan,diff(ring));
%             nE = sum(abs(diffring)>1);
%             disp(nE);
            %diffring(abs(diffring)>1)=nan; % nan out multiple ring crossings
            diffring(diffring<-1)=-1;
            diffring(diffring>1)=1;
            
            for iR=2:nR-1
                k = ring==iR;
                dk = abs(diffring(k))>0;
                pIn(iR,iM,iS) = sum(diffring(k)==-1)./sum(dk);
                pOut(iR,iM,iS) = sum(diffring(k)==1)./sum(dk);
            end
            
            % repeat for random
            [indx0,ring]=find(inRrand==1);
            [~,indx1] = sort(indx0);
            ring = ring(indx1);
            
            diffring = cat(1,nan,diff(ring));
            diffring(diffring>1)=1; % see if this works...
            diffring(diffring<-1)=-1;
            
            for iR=1:nR
                k = ring==iR;
                dk = abs(diffring(k))>0;
                pInRand(iR,iM,iS) = sum(diffring(k)==-1)./sum(dk);
                pOutRand(iR,iM,iS) = sum(diffring(k)==1)./sum(dk);
            end
            
        end
    end
end
radii=radius./11.2;

