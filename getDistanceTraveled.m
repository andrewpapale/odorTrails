% get distance traveled / as the crow flies (straight-line) distance

nM = max(mouse);
nS = max(sess);


Bdistrat = nan(nM,nS);
Ndistrat = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0 & t2s(iM,iS)>0 %#ok<AND2>
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            
            nearX0 = nearX(k);
            nearY0 = nearY(k);
            
            keep = zeros(size(x0));
            frame0 = round(t2s(iM,iS)*50);
            keep(1:frame0)=1;
            notedge = nearX0 > 44.8 & nearX0 < (1280-44.8) & nearY0 > 44.8 & nearY0 < (1024-44.8);
            keep(~notedge)=0;
        
            Bdist = nansum(sqrt(x0(keep==1).^2+y0(keep==1).^2))./11.2;
            Ndist = nansum(sqrt(nx0(keep==1).^2+ny0(keep==1).^2))./11.2;
            
            Bdistrat(iM,iS)=Bdist;%./initDb(iM,iS);
            Ndistrat(iM,iS)=Ndist;%./initD(iM,iS);
            
        end
    end
    disp(iM);
end

zBdr = nan(iM,iS);
zNdr = nan(iM,iS);
for iM=1:nM
    keep = mouse==iM;
    zBdr(iM,:)=nanzscore(Bdistrat(iM,:));
    zNdr(iM,:)=nanzscore(Ndistrat(iM,:));
end
            