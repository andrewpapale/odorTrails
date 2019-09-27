% divide arenas into N partitions
% compute fraction of arena explored / session
% compute number of revisits of each partition
% repeat over different sized partitions

maxX = 55;
maxY = 42;

nM = max(mouse);
nS = max(sess);

parts = (3:100).^2;
nP = length(parts);

% minreturntime = [3,6,9,12,15,18,21,24,27,30];
% nT = length(minreturntime);

pExplored = nan(nP,nM,nS);
% pReturned = nan(nM,nS,nP,nT);
% tReturned = nan(nM,nS,nP,nT);
fExplored = nan(nM,nS);
% fReturned = nan(nM,nS,nP,nT);
for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS;
        kall = sum(k);
        
        if kall > 0
            
            fExplored(iM,iS) = sum(k)./sum(kall);
            
            parfor iN=1:nP % number of partitions
                
                nbinsx = linspace(-55,55,ceil(sqrt(parts(iN))));
                nbinsy = linspace(-42,42,ceil(sqrt(parts(iN))));
                
                [H,~,~,Xct,Yct]=histcounts2(x(k),y(k),nbinsx,nbinsy);
                
                H0 = H(:);
                pExplored(iN,iM,iS) = sum(H0>0)./length(H0);
                
                
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
        end
        disp(iS);
    end
    disp(iM);
end