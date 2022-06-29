%  p(X|X)

% construct tuning curves
nM = max(mouse);
nS = max(sess);
nR = 50;
nTh = 45;
nT = 1000;
rB = quantile(dnT,nR);
tB = quantile(abs(theta),nTh);
timeB = linspace(1,60*50,nT);
iC=1;
TC = nan(1000,nR*nTh);
Q = nan(1000,nT);
for iM=1:nM
    for iS=1:nS
        k = sess==iS & mouse==iM;
        
        if sum(k)>0
            k1 = k & ~edge;
            Q0 = k1 & angcast;
            if any(Q0)
                H = histcn([dnT(Q0),abs(theta(Q0))],rB,tB);
                Occ = histcn([dnT(k1),abs(theta(k1))],rB,tB);
                
                if size(H,1)==nR-1 && size(H,2)==nTh
                    H = cat(1,H,nan(1,nTh));
                end
                if size(Occ,1)==nR-1 && size(Occ,2)==nTh
                    Occ = cat(1,Occ,nan(1,nTh));
                end
                if size(H,2)==nTh-1 && size(H,1)==nR
                    H = cat(2,H,nan(nR,1));
                end
                if size(Occ,2)==nTh-1 && size(Occ,1)==nR
                    Occ = cat(2,Occ,nan(nR,1));
                end
                if size(H,1)==nR-1 && size(H,2)==nTh-1
                    H = cat(1,H,nan(1,nTh-1));
                    H = cat(2,H,nan(nR,1));
                end
                if size(Occ,1)==nR-1 && size(Occ,2)==nTh-1
                    Occ = cat(1,Occ,nan(1,nTh-1));
                    Occ = cat(2,Occ,nan(nR,1));
                end
                Hrt0 = H./(Occ+1E-100);
            else
                Hrt0 = zeros(nR,nTh);
            end
            if any(k1)
                Hq = hist(Q0,timeB);
            else
                Hq = zeros(length(timeB),1);
            end
            TC(iC,1:length(Hrt0(:)))=Hrt0(:);
            Q(iC,:)=Hq/50;
            iC=iC+1;
        end
    end
    disp(iM);
end


shape = size(TC);
nB = prod(shape(2:end));
Px = 1/nB;
nT = size(Q,2);

pxs = nan(nT,nB);

for iB=1:nB
        tempProd = nansum(log(repmat(TC(:,iB),1,nT).^Q));
        tempSum = exp(-0.02*nansum(TC(:,iB)));
        pxs(:,iB) = exp(tempProd)*tempSum*Px;     
        disp(iB);
end

for iT = 1:nT
    pxs(iT,:) = pxs(iT,:)./nansum(pxs(iT,:));
end

pxs = reshape(pxs,[nT,nR,nTh]);   


