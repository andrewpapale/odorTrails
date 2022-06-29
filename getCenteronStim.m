% getCenterOnStim
% 2019-07-26-AndyP
nM = max(mouse);
nS = max(sess);
dt = 3;
tbins = 0.1;
xbins = linspace(-200,200,100);
ybins = linspace(-200,200,80);
PETHt = nan(length(xbins),length(ybins));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kL= mouseL==iM & sessL==iS;
        if sum(k)>0
            k1 = ~edge(k) & bV(k) < 100 & nV(k) < 100 & t2s2(k)==1;
            nx0 = nx(k);
            ny0 = ny(k);
            nV0 = nV(k);
            nx0(k1)=nan;
            ny0(k1)=nan;
            t0 = t(k);
            frame0 = frame(k);
            timOut0 = timOut(kL);
            fraOut0 = fraOut(kL);
            for iT=1:sum(kL)
                zeroframe = interp1(t0,frame0,timOut0(iT),'nearest');
                t00 = unique(interp1(t0,t0,timOut0(iT):0.015:timOut0(iT)+dt,'nearest'));
                f11 = unique(interp1(t0,frame0,timOut0(iT):0.015:timOut0(iT)+dt,'nearest'));
                f11(isnan(f11))=[];
                center = find(f11==zeroframe);
                
                xc0 = nx0(f11)-nx0(fraOut0(iT));
                yc0 = ny0(f11)-ny0(fraOut0(iT));
                nV1 = nV0(f11);
                
                Pt0 = nan(size(PETHt));
                H = histcn([xc0,yc0],xbins,ybins,'AccumData',nV1,'fun',@nanmean);
                if size(H,2)<80
                    H(:,end+1)=nan;
                end
                if size(H,1)<100
                    H(end+1,:)=nan;
                end
                H(H==0)=nan;
                PETHt = cat(3,PETHt,H);
                %Nt = cat(4,Nt,N);
                
            end
        end
        disp(iS);
    end
    disp(iM);
end