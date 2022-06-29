Dline = [];
slope = [];
theta = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            time0 = t(k);
            
            nT = sum(k);
            nP = nT;
            dx = foaw_diff_varTs(x0(1:nP),time0(1:nP), 50, 0.2,0.1)./pixpercm;
            dy = foaw_diff_varTs(y0(1:nP),time0(1:nP), 50, 0.2,0.1)./pixpercm;
            
            Dline0 = nan(nT,1);
            theta0 = nan(nT,1);
            slope0 = nan(nT,1);
            
            for iT=4:nT
                dy0 = (y0(iT)+nanmean(dy(iT-3:iT).*pixpercm))-y0(iT);
                dx0 = (x0(iT)+nanmean(dx(iT-3:iT).*pixpercm))-x0(iT);
                m0 = dy0./dx0;
                slope0(iT)=m0;
                minX = 1;
                maxX = 576;
                result = computeline([x0(iT), y0(iT)],m0, [minX maxX]);
                nP = length(result);
                x1 = result{1}(:,1);
                y1 = result{1}(:,2);
                x2 = result{nP}(:,1);
                y2 = result{nP}(:,2);
                
                Dline0(iT) = ((y2-y1)*nx0(iT)-(x2-x1)*ny0(iT)+x2*y1-y2*x1)./sqrt((y2-y1).^2+(x2-x1).^2);
                dxline = sqrt((nx0(iT)-x0(iT)).^2+(ny0(iT)-x0(iT)).^2);
                theta0(iT) = acosd(dxline/sqrt(dxline.^2+Dline0(iT).^2));
            end
            slope = cat(1,slope,slope0);
            Dline = cat(1,Dline,Dline0./pixpercm);
            theta = cat(1,theta,theta0);
        end
        disp(iS);
    end
end