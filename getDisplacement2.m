% get displacement
M = size(t2s);
dn = nan(M);
db = nan(M);
nbrat = nan(M);
ninitrat = nan(M);
%w0 = nan(M);
sp0 = nan(M);
tnum = nan(M);
Idphi = nan(M);
ndefl = nan(M);
%wp0tr = nan(M);
mdphi = nan(M);
nV1 = nan(M);
m2 = nan(M);
s2 = nan(M);
initD = nan(5,105);
conc2 = nan(M);

nM = 5;
nS = 105;
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS & confn > 0.3 & confb > 0.3 & t2s2==1;
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
            nV0 = nV(k);
            %wp0tr0 = wP(find(k==1,1,'first'));
            for iT=1:ntr0
                m2(iC) = iM;
                s2(iC) = iS;
                dn(iC)=nansum(sqrt((nx0(tr0==iT)).^2+(ny0(tr0==iT)).^2)./5.174);
                db(iC)=nansum(sqrt((x0(tr0==iT)).^2+(y0(tr0==iT)).^2)./5.174);
                nbrat(iC) = dn(iC)./db(iC);
                conc2(iC)=C(iM,iS);
                initD(iM,iS) = sqrt((cx0(iM,iS)-575).^2+(cy0(iM,iS)-15).^2)./5.174; 
                
                ninitrat(iC) = dn(iC)./initD(iM,iS);
                %w0(iC) = Wind(iM,iS);
                %wp0tr(iC) = wp0tr0;
                sp0(iC) = sP(iM,iS);
                nV1(iC) = nanmean(nV0(tr0==iT));
                tnum(iC) = iT;
                Idphi(iC) = nansum(dphi0(tr0==iT));
                mdphi(iC) = nanmean(dphi0(tr0==iT));
                ndefl(iC) = nanmean(abs(Dline0(tr0==iT)));
                iC = iC+1;
            end
        end
    end
    disp(iM);
end
            
            