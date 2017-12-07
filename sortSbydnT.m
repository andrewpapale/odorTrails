% sortSbydnT

nB = 30;
bdnT = quantile(Sdnt,nB);
%bdnT = linspace(0,90,nB);
Ssort = [];
for iB=1:nB
    if iB==1
        k1 = find(Sdnt>0 & Sdnt<=bdnT(iB));
    elseif iB==nB
        k1 = find(Sdnt>bdnT(iB) & Sdnt<=Inf);
    else
        k1 = find(Sdnt>bdnT(iB) & Sdnt<=bdnT(iB+1));
    end
    
    S2 = nanmean(abs(S(k1,:)),1); 
    Ssort = cat(1,Ssort,S2);
end