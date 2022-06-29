minT = 0.1;
nT = length(minT);
nP = ceil(sqrt(nT));
F = figure(1); clf;
for iT=1:nT
    tempAL = nan(15,1);
    tempKP = nan(15,1);
    tempALb = nan(15,1);
    tempKPb = nan(15,1);
    for iR=1:15
        k = ThrV==iR & bait0==0;
        k1 = ThrV==iR & bait0==1;
        tempAL(iR,1)=nansum(t2s(k & ALorKP0==0)>0 & dwellT(k & ALorKP0==0)>minT(iT))./nansum(k(:) & ALorKP0(:)==0);
        tempKP(iR,1)=nansum(t2s(k & ALorKP0==1)>0 & dwellT(k & ALorKP0==1)>minT(iT))./nansum(k(:) & ALorKP0(:)==1);
        tempALb(iR,1)=nansum(t2s(k1 & ALorKP0==0)>0 & dwellT(k1 & ALorKP0==0)>minT(iT))./nansum(k1(:) & ALorKP0(:)==0);
        tempKPb(iR,1)=nansum(t2s(k1 & ALorKP0==1)>0 & dwellT(k1 & ALorKP0==1)>minT(iT))./nansum(k1(:) & ALorKP0(:)==1);
        
    end
    subplot(nP,nP,iT);
    plot(thrVec,tempAL,'linewidth',2); hold on;
    plot(thrVec,tempKP,'--','linewidth',2);
    plot(thrVec,tempALb,'linewidth',2);
    plot(thrVec,tempKPb,'--','linewidth',2);
    if iT==1
        legend('AL unbaited','KP unbaited','AL baited','KP baited');
    end
    title(sprintf('minT = %0.2f',minT(iT)));
    set(gca,'YLim',[0 1]);
end

F = figure(2); clf;
for iT=1:nT
    tempAL0 = nan(15,1);
    tempKP0 = nan(15,1);
    tempAL1 = nan(15,1);
    tempKP1 = nan(15,1);
    tempAL2 = nan(15,1);
    tempKP2 = nan(15,1);
    tempALb = nan(15,1);
    tempKPb = nan(15,1);
    for iR=1:15
        k0 = ThrV==iR & bait0==0 & conc0==0;
        k1 = ThrV==iR & bait0==0 & conc0==1;
        k2 = ThrV==iR & bait0==0 & conc0==2;
        kb = ThrV==iR & bait0==1;
        %k1 = ThrV==iR & bait0==1;
        tempAL0(iR,1)=nansum(t2s(k0 & ALorKP0==0)>0 & dwellT(k0 & ALorKP0==0)>minT(iT))./nansum(k0(:) & ALorKP0(:)==0);
        tempKP0(iR,1)=nansum(t2s(k0 & ALorKP0==1)>0 & dwellT(k0 & ALorKP0==1)>minT(iT))./nansum(k0(:) & ALorKP0(:)==1);
        tempAL1(iR,1)=nansum(t2s(k1 & ALorKP0==0)>0 & dwellT(k1 & ALorKP0==0)>minT(iT))./nansum(k1(:) & ALorKP0(:)==0);
        tempKP1(iR,1)=nansum(t2s(k1 & ALorKP0==1)>0 & dwellT(k1 & ALorKP0==1)>minT(iT))./nansum(k1(:) & ALorKP0(:)==1);
        tempAL2(iR,1)=nansum(t2s(k2 & ALorKP0==0)>0 & dwellT(k2 & ALorKP0==0)>minT(iT))./nansum(k2(:) & ALorKP0(:)==0);
        tempKP2(iR,1)=nansum(t2s(k2 & ALorKP0==1)>0 & dwellT(k2 & ALorKP0==1)>minT(iT))./nansum(k2(:) & ALorKP0(:)==1);
        tempALb(iR,1)=nansum(t2s(kb & ALorKP0==0)>0 & dwellT(kb & ALorKP0==0)>minT(iT))./nansum(kb(:) & ALorKP0(:)==0);
        tempKPb(iR,1)=nansum(t2s(kb & ALorKP0==1)>0 & dwellT(kb & ALorKP0==1)>minT(iT))./nansum(kb(:) & ALorKP0(:)==1);
        
    end
    subplot(nP,nP,iT);
    plot(thrVec,tempAL0,'linewidth',2); hold on;
    plot(thrVec,tempKP0,'--','linewidth',2);
    plot(thrVec,tempAL1,'linewidth',2);
    plot(thrVec,tempKP1,'--','linewidth',2);
    plot(thrVec,tempAL2,'linewidth',2);
    plot(thrVec,tempKP2,'--','linewidth',2);
    plot(thrVec,tempALb,'k-','linewidth',2);
    plot(thrVec,tempKPb,'k--','linewidth',2);
    if iT==1
        legend('AL 0.1%','KP 0.1%','AL 1.0%','KP 1.0%','AL 2.0%','KP 2.0%','AL baited','KP baited');
    end
    title(sprintf('minT = %0.2f',minT(iT)));
    set(gca,'YLim',[0 1]);
end

F = figure(3); clf;
iC=1;
minT = 0.5;
for iN=0:2
    for iM=1:nM
        tempAL = nan(15,1);
        tempKP = nan(15,1);
        tempALb = nan(15,1);
        tempKPb = nan(15,1);
        for iR=1:15
            k = ThrV==iR & bait0==0 & mouse0==iM & conc0==iN;
            k1 = ThrV==iR & bait0==1 & mouse0==iM & conc0==iN;
            tempAL(iR,1)=nansum(t2s(k & ALorKP0==0)>0 & dwellT(k & ALorKP0==0)>minT(iT))./nansum(k(:) & ALorKP0(:)==0);
            tempKP(iR,1)=nansum(t2s(k & ALorKP0==1)>0 & dwellT(k & ALorKP0==1)>minT(iT))./nansum(k(:) & ALorKP0(:)==1);
            tempALb(iR,1)=nansum(t2s(k1 & ALorKP0==0)>0 & dwellT(k1 & ALorKP0==0)>minT(iT))./nansum(k1(:) & ALorKP0(:)==0);
            tempKPb(iR,1)=nansum(t2s(k1 & ALorKP0==1)>0 & dwellT(k1 & ALorKP0==1)>minT(iT))./nansum(k1(:) & ALorKP0(:)==1);
            
        end
        subplot(3,nM,iC);
        plot(thrVec,tempAL,'linewidth',2); hold on;
        plot(thrVec,tempKP,'--','linewidth',2);
        plot(thrVec,tempALb,'linewidth',2);
        plot(thrVec,tempKPb,'--','linewidth',2);
        if iM==1 && iN==1
            legend('AL unbaited','KP unbaited','AL baited','KP baited');
        end
        %title(sprintf('minT = %0.2f',minT(iT)));
        set(gca,'YLim',[0 1]);
        iC=iC+1;
    end
end