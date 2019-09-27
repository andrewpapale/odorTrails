%shuffle spots
nSh = 100;
spotThr = 5; % cm
spotThr = spotThr*11.2; % cm->pixels


nM = max(mouse);
nS = max(sess);
xT = [];
yT = [];

for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS & ~removeint2s & ~edgespot;
        
        if sum(k)>0
            
            kT = mouseT==iM & sessT==iS;
            
            % get all spots
            xT = cat(1,xT,nanmedian(xT1(kT)));
            yT = cat(1,yT,nanmedian(yT1(kT)));
            
        end
    end
    disp(iM);
end

%% shuffle spots
rng('shuffle');
D = pdist2([xT,yT],[xT,yT],'Euclidean');
D(D<spotThr)=nan;

Hori = [];
Hdphi = [];
Hocc_not = [];
Hocc_yes = [];
Hps = [];
Hnv = [];
sminD = nan(nM,nS,nSh);
sht2s = nan(nM,nS,nSh);
sinitD = nan(nM,nS,nSh);
sDnid = nan(nM,nS,nSh);
origsucc = nan(nM,nS,nSh);
xS = nan(nM,nS,nSh);
yS = nan(nM,nS,nSh);

for iSh = 1:nSh
    xT2 = nan(size(xT));
    yT2 = nan(size(yT));
    for iT=1:length(xT)
        row = D(iT,:);
        k0 = find(~isnan(row));
        if ~isempty(k0)
            six = randi(length(k0));
            xT2(iT) = xT(k0(six));
            yT2(iT) = yT(k0(six));
            D(six,:) = nan;
        else
            xT2(iT) = nan;
            yT2(iT) = nan;
        end
    end
    
    sdnt = nan(size(x));
    orients = nan(size(x));
    st2s = zeros(size(x));
    st2s1 = nan(size(cx0));
    Dnid = nan(size(cx0));
    Dnd = nan(size(cx0));
    sspotfound = zeros(size(x));
    sspotfound0trunc = zeros(size(x));
    lastt2s = 1500;
    for iM=1:nM
        for iS=1:nS
            
            k = mouse==iM & sess==iS;
            
            if sum(k)>0 & ~isnan(xT2(iT)) %#ok<AND2>
                
                x0 = x(k);
                y0 = y(k);
                nx0 = nx(k);
                ny0 = ny(k);
                frame0 = frame(k);
                nV0 = nV(k);
                V0 = V(k);
                edge0 = edge(k);
                dphi0 = log10dphi(k);
                orient0 = orient(k);
                
                % calculate distance from trail
                dnT0 = sqrt((nx0-xT2(iT)).^2+(ny0-yT2(iT)).^2);
                sdnt(k) = dnT0;
                sminD(iM,iS,iSh) = nanmin(dnT0);
                foundspot = find(dnT0 < threshold,1,'first');
                
                BSx = xT2(iT)-x0;
                BSy = yT2(iT)-y0;
                BNx = nx0-x0;
                BNy = ny0-y0;
                
                orient0s = wrapTo180(atan2d(BNy,BNx)-atan2d(BSy,BSx));
                orients = cat(1,orients,orient0s);
                
                st2s0 = zeros(size(dnT0));
                
                
                if ~isempty(foundspot)
                    lastt2s = foundspot;
                    sspotfound(k)=ones(sum(k),1);
                    k0 = find(k==1,1,'first'):find(k==1,1,'first')+foundspot;
                    st2s(k0)=ones(length(k0),1);
                    st2s1(iM,iS) = foundspot/50;
                    st2s0(1:foundspot)=1;
                    k1 = st2s0 & ~edge0 & nV0 < 100 & V0 < 100;
                else
                    k0 = find(k==1,1,'first'):find(k==1,1,'last');
                    sspotfound0trunc(k0) = ones(length(k0),1);
                    st2s1(iM,iS)=-20;
                    k1 = ~edge0 & nV0 < 100 & V0 < 100;
                end
                
                origsucc(iM,iS,iSh)=spotfound(iM,iS);
                
                Dnd(iM,iS) = nansum(sqrt(diff(nx0(k1)).^2+diff(ny0(k1)).^2)./11.2);% displacement of nose
                k2 = find(~isnan(dnT0) & frame0 <= 50,1,'first');
                sht2s(iM,iS,iSh) = st2s1(iM,iS);
                if ~isempty(k2)
                    sinitD(iM,iS,iSh) = dnT0(k2);
                    sDnid(iM,iS,iSh)= Dnd(iM,iS)./sinitD(iM,iS,iSh); % displacement of nose/initial distance to spot
                end
 
                bins = linspace(1.5,100,30);
                k3 = ~edge(k) & ~isstopped(k) & (t2s0(k)==1 | spotfound0trunc(k)==1) & nV0 < 100 & V0 < 100;
                Hdphi0 = histcn(dnT0(k3),bins,'AccumData',dphi0(k3),'fun',@nanmean);
                Hdphi(1:29,iM,iS,iSh) = Hdphi0(1:29,1,1,1);
                Hnv0 = histcn(dnT0(k3),bins,'AccumData',nV0(k3),'fun',@nanmean);
                Hnv(1:29,iM,iS,iSh)=Hnv0(1:29,1,1,1);
                
                xS(iM,iS,iSh) = xT2(iT);
                yS(iM,iS,iSh) = yT2(iT);
                
            end
        end
        disp(iM);
    end
    
    
    
    
    % occupancy
    radius = linspace(1.5,100,100);
    nR = length(radius);
    area = nan(nR,1);
    r0 = cat(2,0,radius);
    for iR=2:nR
        area(iR) = pi*(r0(iR).^2-r0(iR-1).^2);
    end
    nM = 5; nS = 197;
    R = nan(nR,nM,nS);
    for iM=1:nM
        for iS=1:nS
            k = mouse==iM & sess==iS & nV < 100 & V < 100;
            kT = mouseT==iM & sessT==iS;
            if sum(k)>0
                skipflag = 0;
                if st2s1(iM,iS)>mint2s
                    k0 = k & st2s==1;% & ~isstopped;
                elseif st2s1(iM,iS)==-20
                    k0 = k;
                elseif st2s1(iM,iS)==-40
                    skipflag = 1;
                end
                if ~skipflag
                    temp = histcounts(sdnt(k0),radius);
                    R(1:length(temp),iM,iS)=temp;%./nansum(temp,2);
                end
            end
        end
        disp(iM);
    end
    radii = radius;
    
    minB = 1;
    %area = 1;
    k = st2s1>mint2s;
    pRs = R(:,k);
    N = nansum(pRs,2);
    pRs(N<minB,:)=nan;
    N0 = nansum(pRs,1);
    N0(N0<minB)=nan;
    mpRs = nanmean(pRs./(repmat(N0,[size(pRs,1),1]).*area),2);
    dpRs = nanstderr(pRs./(repmat(N0,[size(pRs,1),1]).*area),[],2);
    k = st2s1==-20;
    pRf = R(:,k);
    N = nansum(pRf,2);
    pRf(N<minB,:)=nan;
    N0 = nansum(pRf,1);
    N0(N0<minB)=nan;
    mpRf = nanmean(pRf./(repmat(N0,[size(pRf,1),1]).*area),2);
    
    Hocc_not = cat(2,Hocc_not,mpRf);
    Hocc_yes = cat(2,Hocc_yes,mpRs);
    
    
    hit = nan(length(radii),nM);
    miss = nan(length(radii),nM);
    kh = st2s1>mint2s;
    km = st2s1==-20;
    for iM=1:nM
        kh0 = kh & mouse0==iM;
        km0 = km & mouse0==iM;
        for iR=1:length(radii)
            sRh = nansum(R(iR,kh));
            sRm = nansum(R(iR,km));
            hit(iR,iM) = nansum(R(iR,kh0));
            miss(iR,iM) = nansum(R(iR,km0));
            
        end
    end
    pSR = (hit./(hit+miss));
    
    Hps = cat(2,Hps,nanmean(pSR,2));
    
    k = ~edge & ~isstopped & t2s0==1 & nV < 100 & V < 100;
    Hori0 = histogram2(sdnt(k),abs(orients(k)),linspace(1.5,100,50),linspace(0,180,45),'facecolor','flat');
    Hori(:,:,iSh) = Hori0;
    
    
    fprintf('shuffle %d/%d \n',iSh,nSh);
end
