% TimeInRing
% 2017-08-07 AndyP
% get the time in each ring around spots

nR = 100;
dR = 5.174; % cm/pixel

nM = max(mouse);
nS = max(sess);

R = nan(nR,nM,nS);
Rrat = nan(nR,nM,nS);
randR = nan(nR,nM,nS);
radii0 = nan(nR,nM,nS);
% sess0 = nan(nR,nM,nS);
% trial0 = nan(nR,nM,nS);
% conc0 = nan(nR,nM,nS);
% mouse0 = nan(nR,nM,nS);
% ALorKP0 = nan(nR,nM,nS);
% bait0 = nan(nR,nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            xT = cx0(iM,iS);
            yT = cy0(iM,iS);
            
            x0 = nx(k);
            y0 = ny(k);
            edge0 = edge(k);
            
            x0 = x0(~edge0);
            y0 = y0(~edge0);
            
            
            
            % rand data for shuffle control
            xR = nanmin(x)+(nanmax(x)-nanmin(x)).*rand(length(x0),1);
            yR = nanmin(y)+(nanmax(y)-nanmin(y)).*rand(length(y0),1);
            
            % create circles around spot
            
            
            dist = sqrt((x0-xT).^2+(y0-yT).^2);
            randdist = sqrt((xR-xT).^2+(yR-yT).^2);
            radius = dR.*(1:nR);
            for iR=1:nR
                if iR>1
                    inR = dist <= radius(iR) & dist > radius(iR-1);
                    inRrand = randdist <= radius(iR) & randdist > radius(iR-1);
                    area = pi*(radius(iR).^2-radius(iR-1).^2);
                elseif iR==1
                    inR = dist <= radius(iR);
                    inRrand = randdist <= radius(iR);
                    area = pi*radius(iR).^2;
                end
                
                %area=1;
                R(iR,iM,iS) = nansum(inR);%./(area);
                Rrat(iR,iM,iS) = nansum(inR)./(area*length(x0));
                
                randR(iR,iM,iS) = nansum(inRrand)./(area*length(x0));
            end
            
            
            radii0(:,iM,iS)=radius./5.174;
            %             conc0(:,iM,iS)=repmat(mode(conc(k)),[nR,1]);
            %             sess0(:,iM,iS)=repmat(mode(sess(k)),[nR,1]);
            %             mouse0(:,iM,iS)=repmat(mode(mouse(k)),[nR,1]);
            %             ALorKP0(:,iM,iS)=repmat(mode(ALorKP(k)),[nR,1]);
            %             trial0(:,iM,iS)=repmat(mode(trial(k)),[nR,1]);
            %             bait0(:,iM,iS)=repmat(mode(bait(k)),[nR,1]);
        end
    end
    disp(iM);
end
radii=radii0(:,1)';



% minB = 1;
% area = 1;
% pRs = R(:,cond==1);
% N = nansum(pRs,2);
% pRs(N<minB,:)=nan;
% N0 = nansum(pRs,1);
% N0(N0<minB)=nan;
% mpRs = nanmean(pRs./(repmat(N0,[size(pRs,1),1]).*area),2);
% dpRs = nanstderr(pRs./(repmat(N0,[size(pRs,1),1]).*area),[],2);
% pRf = R(:,cond==0);
% N = nansum(pRs,2);
% pRf(N<minB,:)=nan;
% N0 = nansum(pRf,1);
% N0(N0<minB)=nan;
% mpRf = nanmean(pRf./(repmat(N0,[size(pRf,1),1]).*area),2);
% dpRf = nanstderr(pRf./(repmat(N0,[size(pRf,1),1]).*area),[],2);
% 
% F = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial'); 
% clf;
% 
% lineProps.col = {'r'};
% mseb(radii(1:end),mpRs(1:end)',dpRs(1:end)',lineProps);
% %mseb(radii,log10(mpRs'),dpRs'./mpRs'*log(10),lineProps);
% lineProps.col = {'b'};
% % mseb(radii,log10(mpRf'),dpRf'./mpRf'*log(10),lineProps);
% mseb(radii(1:end),mpRf(1:end)',dpRf(1:end)',lineProps);
% set(gca,'fontsize',36);
% xlabel('Distance From Spot (cm)','fontsize',24,'FontName','Arial');
% ylabel('Percent Time (%)','fontsize',24,'FontName','Arial');






% L = line([1.5 1.5],[0 0.05]);
% set(L,'LineStyle','--','color','k','LineWidth',2);
% % 
% L1 = legend('Spot Found','CSM','CRW Control');
% set(L1,'fontsize',24,'FontName','Arial');
% 
% %set(gca,'YLim',[0 0.04],'YTick',[0 0.005 0.01, 0.015, 0.02, 0.025 0.03 0.035 0.04],'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'});
% set(gca,'FontName','Arial');




