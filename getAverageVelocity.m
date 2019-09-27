% get average velocity for each session

nM = max(mouse);
nS = max(sess);

nV1 = nan(nM,nS);
zdp1 = nan(nM,nS);
%edge = ~(nearX>50 & nearY>50 & nearX<1280-50 & nearY<1024-50);
%t2s0 = 30*50;
for iM=1:nM
    for iS=1:nS
        keep = mouse==iM & sess==iS;
        
        if sum(keep)>0
            
            nV0 = nV(keep & ~edge & dnT<100);
            znidphi3 = znidphi2(keep & ~edge);
%             if t2s(iM,iS)>0
%                 t2s0 = round(t2s(iM,iS)*50);
%             else
%             end
%                 nV0 = nV0(1:t2s0);
                nV1(iM,iS) = nanmean(nV0);
                zdp1(iM,iS) = nanmean(znidphi3);
            
        end
    end
    disp(iM);
end

znV1 = nan(nM,nS);
for iM=1:nM
    znV1(iM,:)=nanzscore(nV1(iM,:));
end