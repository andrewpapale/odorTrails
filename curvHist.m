% curvature histogram


binsx = linspace(-100,500,100);
binsy = linspace(-400,400,150);


H = nan(length(binsx),length(binsy));
N = nan(length(binsx),length(binsy));

k = ~edge & bV < 100 & nV < 100 & ~remove & t2s1;

for ibx = 1:length(binsx)-1
    for iby = 1:length(binsy)-1
        k0 =k & (binsx(ibx) <= nxcr < binsx(ibx+1)) & (binsx(iby) <= nycr < binsy(iby+1));
        H(ibx,iby) = nanmean(log10dphi(k0));
        N(ibx,iby) = sum(k0);
    end
end