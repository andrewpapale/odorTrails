% get displacement

dn = nan(size(t2s));
db = nan(size(t2s));
nbrat = nan(size(t2s));
ninitrat = nan(size(t2s));
w0 = nan(size(t2s));
sp0 = nan(size(t2s));
tnum = nan(size(t2s));
Idphi = nan(size(t2s));
ndefl = nan(size(t2s));
wp0tr = nan(size(t2s));
mdphi = nan(size(t2s));
typeS = nan(size(t2s));
ampS = nan(size(t2s));
initDS = nan(size(t2s));
concS = nan(size(t2s));
LapNo = nan(size(t2s));
nM = 3;
nS = 40;
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS & ~edge & t2s2==1;
        kS = mouseT==iM & sessT==iS;
        kL = mouseL==iM & sessL==iS;
        if sum(k)>0
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            dphi0 = log10dphi(k);
            Dline0 = Dline(k);
            dnT0 = dnT(k);
            tr0 = t2s0(k);
            ntr0 = max(tr0);
            wp0tr0 = wP(find(k==1,1,'first'));
            tff0 = tff(kS);
            t2s3 = t2s(kS);
            ampOut0 = ampOut(kL);
            timeOut0 = timOut(kL);
            for iT=1:ntr0
                ampS(iC) = mode(ampOut0(timeOut0 > tff0(iT) & timeOut0 < t2s3(iT)));
                dn(iC)=nansum(sqrt((nx0(tr0==iT)).^2+(ny0(tr0==iT)).^2)./5.174);
                db(iC)=nansum(sqrt((x0(tr0==iT)).^2+(y0(tr0==iT)).^2)./5.174);
                nbrat(iC) = dn(iC)./db(iC);
                LapNo(iC) = iT;
                initDS(iC) = initD(iM,iS);
                ninitrat(iC) = dn(iC)./initD(iM,iS);
                w0(iC) = Wind(iM,iS);
                wp0tr(iC) = wp0tr0;
                sp0(iC) = sP(iM,iS);
                tnum(iC) = iT;
                Idphi(iC) = nansum(dphi0(tr0==iT & dnT0 < 50));
                mdphi(iC) = nanmean(dphi0(tr0==iT & dnT0 < 50));
                ndefl(iC) = nanmean(abs(Dline0(tr0==iT & dnT0 < 50)));
                typeS(iC) = cond(iM,iS);
                concS(iC) = C(iM,iS);
                iC = iC+1;

            end
        end
    end
    disp(iM);
end
            
            