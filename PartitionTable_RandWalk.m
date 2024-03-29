% divide arenas into N partitions
% compute fraction of arena explored / session
% compute number of revisits of each partition
% repeat over different sized partitions

maxX = 1280;
maxY = 1020;

nM = max(mouse);
nS = max(sess);

parts = (3:100).^2;
nP = length(parts);

minreturntime = [3,6,9,12,15,18,21,24,27,30];
nT = length(minreturntime);

pExp_RW = nan(nM,nS,nP);
pRet_RW = nan(nM,nS,nP,nT);
tRet_RW = nan(nM,nS,nP,nT);
fExp_RW = nan(nM,nS);
fRet_RW = nan(nM,nS,nP,nT);
for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS & ~edge & V<100;
        kall = sum(mouse==iM & sess==iS);
        
        fExp_RW(iM,iS) = sum(k)./sum(kall);
        
        for iN=1:nP % number of partitions
            
            nbinsx = linspace(44.8,maxX-44.8,ceil(sqrt(parts(iN))));
            nbinsy = linspace(44.8,maxY-44.8,ceil(sqrt(parts(iN))));
            
            [H,~,~,Xct,Yct]=histcounts2(x(k),y(k),nbinsx,nbinsy);
            
            H0 = H(:);
            pExp_RW(iM,iS,iN) = sum(H0>0)./length(H0);
            
            
            % find returns to partition
%             [i,j]=ind2sub(size(H),1:length(H0));
%             for iT=1:nT
%                 sumreturn = 0;
%                 sreturns = 0;
%                 treturn = [];
%                 for iP=1:length(H0)
%                     if H0(iP)>1 % potential return
%                         kxy = find(Xct==i(iP) & Yct==j(iP));
%                         dkxy = cat(1,0,diff(kxy));
%                         anygood = dkxy>=minreturntime(iT)*50;
%                         if any(anygood)
%                             sreturns = sum(anygood) + sreturns;
%                             sumreturn = sumreturn + 1;
%                             treturn = cat(1,treturn, dkxy(anygood)/50);
%                         end
%                     end
%                 end
%                 
%                 pReturned(iM,iS,iN,iT) = sumreturn./length(H0); % percent of partitions returned to
%                 tReturned(iM,iS,iN,iT) = nanmedian(treturn); % median time of returns
%                 fReturned(iM,iS,iN,iT) = sreturns./sum(H0); % percent of time spent re-exploring
%             end
        end
        disp(iS);
    end
    disp(iM);
end