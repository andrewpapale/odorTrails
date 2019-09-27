% RingTransitionMatrix
% 2017-08-08 AndyP
% get transition matrix for found/not found trials

doTest = false;
nR = 31;
dR = 11.2; % cm/pixel

nM = max(mouse);
nS = max(sess);

iC=1;
iD=1;
iE=1;
foundR = [];
notR = [];
ktime = round(nanmedian(t2s(t2s>0))*50);
Cfound = [];
Cnot = [];
Crand = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = nx(k);
            y0 = ny(k);
            
            nearX0 = nearX(k);
            nearY0 = nearY(k);
            
            edge = ~(nearX0>50 & nearY0>50 & nearX0<1280-50 & nearY0<1024-50);
            
            %disp(size(x0));
            
            kyes = t2s(iM,iS)>0;
            kno = t2s(iM,iS)==-20;
            
            if kyes | kno
                
                if kyes
                    ktime = round(t2s(iM,iS)*50)+500;
                else % do not update ktime
                end
                
                %disp(min(length(x0),ktime));
                
                x0 = x0(1:min(length(x0),ktime));
                y0 = y0(1:min(length(y0),ktime));
                edge = edge(1:min(length(x0),ktime));
                
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
                        inR = dist <= radius(iR) & dist > radius(iR-1) & ~edge;
                        inRrand = randdist <= radius(iR) & randdist > radius(iR-1);
                        area = pi*(radius(iR).^2-radius(iR-1).^2);
                    else
                        inR = dist <= radius(iR) & ~edge;
                        inRrand = randdist <= radius(iR);
                        area = pi*radius(iR).^2;
                    end
                    
                    %area=1;
                    
                    keep = ~isnan(inR);
                    
                    if kyes
                        if sum(keep)>0
                            foundR(iR,iC) = nansum(inR(keep))./(50*area);
                        else
                            foundR(iR,iC)=nan;
                        end
                    elseif kno
                        if sum(keep)>0
                            notR(iR,iD)=nansum(inR(keep))./(50*area);
                        else
                            notR(iR,iD)=nan;
                        end
                        keep = ~isnan(inRrand);
                        if sum(keep)>0
                            randR(iR,iE) = nansum(inRrand(keep))./(50*area);
                        else
                            randR(iR,iE)=nan;
                        end
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
                iE=iE+1;
            end
        end
        
    end
    disp(iM);
end
radii=radius./(11.2);

for iR=1:nR
    for iQ=1:nR
                Cfound(iR,iQ)=corr2(foundR(iR,:),foundR(iQ,:));
                %Cnot(iR,iQ) = corr2(notR(iR,:),notR(iQ,:));
                %Crand(iR,iQ)= corr2(randR(iR,:),randR(iQ,:));
%         Cfound(iR,iQ) = mi(foundR(iR,:)',foundR(iQ,:)');
%         Cnot(iR,iQ) = mi(notR(iR,:)',notR(iQ,:)');
%         Crand(iR,iQ)=mi(randR(iR,:)',randR(iQ,:)');
    end
end

F = figure(1); clf;
subplot(1,3,1);
pcolor(radii,radii,Cfound); shading flat;
colorbar;
%caxis([-0.5 1]);
caxis([0 4]);
title('Spot found','fontsize',21);
xlabel('radius (cm)','fontsize',21);
ylabel('radius (cm)','fontsize',21);
set(gca,'fontsize',18);
% subplot(1,3,2);
% pcolor(radii,radii,Cnot); shading flat;
% %caxis([-0.5 1]);
% caxis([0 4]);
% colorbar;
% title('Spot Missed','fontsize',21);
% xlabel('radius (cm)','fontsize',21);
% set(gca,'fontsize',18);
% subplot(1,3,3);
% pcolor(radii,radii,Crand); shading flat;
% %caxis([-0.5 0.5]);
% caxis([0 4]);
% colorbar;
% title('Random','fontsize',21);
% xlabel('radius (cm)','fontsize',21);
% set(gca,'fontsize',18);

F = figure(2); clf;
plot(radii,nanmean(Cfound,2),'linewidth',2);
hold on;
plot(radii,nanmean(Cnot,2),'linewidth',2);
plot(radii,nanmean(Crand,2),'linewidth',2);
legend('found','missed','random');
xlabel('radius (cm)','fontsize',21);
ylabel('mean MI','fontsize',21);
set(gca,'fontsize',18);