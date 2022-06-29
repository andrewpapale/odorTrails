% TimeInRing3
% 2017-08-07 AndyP
% get the time in each ring around spots

mint2s = 0;
ktime = 1500; % frames, initial condition
radius = linspace(2,100,100);
nR = length(radius);
nM = 3;
nS = 40;
area = nan(nR,1);
for iR=1:nR
    if iR==1
        area(iR)=pi*radius(iR).^2;
    else
        area(iR)=pi*(radius(iR).^2-radius(iR-1).^2);
    end
end
R = nan(nR,nM,nS,10);
AmpR = nan(nM,nS,10);
TypeR = nan(nM,nS,10);
for iM=1:nM
    for iS=1:nS
        nT = sT(iM,iS);
        if sum(k)>0 && ~isnan(nT)
            kL = mouseL==iM & sessL==iS;
            fkL = find(kL);
            ampOut0 = ampOut(kL);
            timeOut0 = timOut(kL);
            t2s1 = t2s(mouseT==iM & sessT==iS);
            tff0 = tff(mouseT==iM & sessT==iS);
            lapt = zeros(size(x));
            dnT0 = dnT(k);
            for iT=1:nT
                k0 = ~edge(k) & mL(k)<25 & nV(k) < 100 & bV(k) < 100 & t2s0(k)==iT;
                temp = histcounts(dnT0(k0),linspace(2,nR,nR));
                R(1:length(temp),iM,iS,iT)=temp;%./nansum(temp,2);
                AmpR(iM,iS,iT)=mode(Amp(k0));
                TypeR(iM,iS,iT)=mode(Type(k0));
            end
        end
    end
    disp(iM);
end
radii = radius;


area = 1;
clf
mint2s = 1;
minB = 1;
pRs = R(:,AmpR==5);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
N0 = nansum(pRs,1);
N0(N0<minB)=nan;
mpRs = nanmean(pRs./(repmat(N0,[size(pRs,1),1]).*area),2);
dpRs = nanstderr(pRs./(repmat(N0,[size(pRs,1),1]).*area),[],2);

lineProps.col = {'r'};
mseb(radii,mpRs',dpRs',lineProps);

pRs = R(:,AmpR==0 & TypeR==1);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
N0 = nansum(pRs,1);
N0(N0<minB)=nan;
mpRs = nanmean(pRs./(repmat(N0,[size(pRs,1),1]).*area),2);
dpRs = nanstderr(pRs./(repmat(N0,[size(pRs,1),1]).*area),[],2);

lineProps.col = {'b'};
mseb(radii,mpRs',dpRs',lineProps);

pRs = R(:,AmpR==0 & TypeR==0);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
N0 = nansum(pRs,1);
N0(N0<minB)=nan;
mpRs = nanmean(pRs./(repmat(N0,[size(pRs,1),1]).*area),2);
dpRs = nanstderr(pRs./(repmat(N0,[size(pRs,1),1]).*area),[],2);

lineProps.col = {'k'};
mseb(radii,mpRs',dpRs',lineProps);
