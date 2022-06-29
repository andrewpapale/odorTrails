% TimeInRing_Y
% 2017-08-07 AndyP
% get the time in each ring around center of Y

nR = 100;
dR = 11.2; % cm/pixel
nM = max(mouse);
nS = max(sess);
minB = 30;

radius = dR.*(1:nR);
area = nan(nR,1);
for iR=1:nR
    if iR==1
        area(iR)=pi*radius(iR).^2;
    else
        area(iR)=pi*(radius(iR).^2-radius(iR-1).^2);
    end
end


R = nan(nR,nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            temp = histcounts(Y.dnTc(k),linspace(0,nR,nR));
            R(1:length(temp),iM,iS)=temp;%./nansum(temp,2);
        end
    end
    disp(iM);
end
radii = radius./11.2;

pRs = R(:,:);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
pRs = pRs'./(nansum(pRs(:)).*area');
mpRs = nanmean(pRs,1);
dpRs = nanstderr(pRs,[],1);

F = figure(1); clf;
lineProps.col = {'r'};
mseb(radii,log10(mpRs'),1./(log10(dpRs)*log(10)),lineProps);
xlabel('Distance From Spot (cm)','fontsize',24);
ylabel('Fractional Occupancy / Area (cm^{-2})','fontsize',24);
set(gca,'fontsize',21);
