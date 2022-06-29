
k = V < 1000 & nV < 1000 & confb > 0.3 & confn > 0.5 & ~(Top & Right) & ~(Top & Left) & ~(Bottom & Right) & ~(Bottom & Left) & ~(~Bottom & ~Top & ~Left & ~Right & ~center) & ~isinf(log10dphiB) & sqrt((nx-x).^2+(ny-y).^2) < 80;

cond = nan(22,1); 
for iS=1:22
    if strcmp(InfoMatrix{iS,3},'A') 
        cond(iS)=0;
    elseif strcmp(InfoMatrix{iS,3},'C') 
        cond(iS)=1; 
    end 
end
Vs = nan(5,22);
for iS=1:22
    k0 = sess==iS & k;
    
%     Vs0(1) = nansum(center(k0));
%     Vs0(2) = nansum(Bottom(k0));
%     Vs0(3) = nansum(Top(k0));
%     Vs0(4) = nansum(Left(k0));
%     Vs0(5) = nansum(Right(k0));
%     %
%     Vs(:,iS)=Vs0./sum(Vs0,2);
    
%     Vs(1,iS) = nanmean(nV(k0 & center)./V(k0 & center));
%     Vs(2,iS) = nanmean(nV(k0 & Bottom)./V(k0 & Bottom));
%     Vs(3,iS) = nanmean(nV(k0 & Top)./V(k0 & Top));
%     Vs(4,iS) = nanmean(nV(k0 & Left)./V(k0 & Left));
%     Vs(5,iS) = nanmean(nV(k0 & Right)./V(k0 & Right));
%     

      Vs(1,iS) = nanmean(log10dphiB(k0 & center));
      Vs(2,iS) = nanmean(log10dphiB(k0 & Bottom));
      Vs(3,iS) = nanmean(log10dphiB(k0 & Top));
      Vs(4,iS) = nanmean(log10dphiB(k0 & Left));
      Vs(5,iS) = nanmean(log10dphiB(k0 & Right));

%     Vs(1,iS) = nanmean(V(k0 & center));
%     Vs(2,iS) = nanmean(V(k0 & Bottom));
%     Vs(3,iS) = nanmean(V(k0 & Top));
%     Vs(4,iS) = nanmean(V(k0 & Left));
%     Vs(5,iS) = nanmean(V(k0 & Right));
%     
   
end

 Vz = nan(5,2); 
 for iZ=1:5
     for iC=0:1 
         Vz(iZ,iC+1) = nanmean(Vs(iZ,cond==iC)); 
         dVz(iZ,iC+1)=nanstderr(Vs(iZ,cond==iC)); 
     end
 end