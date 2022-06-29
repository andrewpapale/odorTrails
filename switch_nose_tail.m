function [tnx,tny,ttx,tty] = switch_nose_tail(x0,y0,time0,nx0,ny0,tx0,ty0)
% flip nose-tail

doTest = false;
if doTest
    figure(1); clf;
end
% dt = cat(1,nan,sqrt(diff(tx0).^2+diff(ty0).^2));
% dn = cat(1,nan,sqrt(diff(nx0).^2+diff(ny0).^2));

% find direction body is traveling
% find minimum point (tx0 or nx0) from vector
dx = foaw_diff_varTs(x0, time0, 50, 0.2,0.1);
dy = foaw_diff_varTs(y0, time0, 50, 0.2,0.1);

nT = length(x0);
tnx = nan(nT,1);
tny = nan(nT,1);
ttx = nan(nT,1);
tty = nan(nT,1);

if doTest
    plot(x0,y0,'.','markersize',1,'color',[0.5 0.5 0.5]);
    hold on;
end

V0 = nan(nT,1);
for iT=4:nT
    
    if ~isnan(x0(iT))
        dy0 = (y0(iT)+nanmean(dy(iT-3:iT)))-y0(iT);
        dx0 = (x0(iT)+nanmean(dx(iT-3:iT)))-x0(iT);
        V = dx0.^2+dy0.^2;
        V0(iT) = V;
        if abs(V) >0.001 && abs(dx0) > 0
            
            m0 = dy0./dx0;
            
            switch sign(dx0)
                case -1
                    minX = x0(iT)-20;
                    maxX = x0(iT);
                case 1
                    minX = x0(iT);
                    maxX = x0(iT)+20;
            end
            
            result = computeline([x0(iT), y0(iT)],m0, [minX maxX]);
            nP = length(result);
            x1 = nan(nP,1);
            y1 = nan(nP,1);
            for iP=1:nP
                x1(iP) = result{iP}(:,1);
                y1(iP) = result{iP}(:,2);
            end
            
            dbnx = nanmin((nx0(iT)-x1).^2);
            dbny = nanmin((ny0(iT)-y1).^2);
            dbn = nanmin(sqrt(dbnx+dbny));
            dbtx = nanmin((tx0(iT)-x1).^2);
            dbty = nanmin((ty0(iT)-y1).^2);
            dbt = nanmin(sqrt(dbtx+dbty));
            
            if dbn>dbt
                tnx(iT) = tx0(iT);
                tny(iT) = ty0(iT);
                ttx(iT) = nx0(iT);
                tty(iT) = ny0(iT);
            else
                tnx(iT) = nx0(iT);
                tny(iT) = ny0(iT);
                ttx(iT) = tx0(iT);
                tty(iT) = ty0(iT);
            end
            
            if doTest
                P7 = plot(x1,y1,'r.-');
                if iT>5
                    P1 = plot(x0(iT-5:iT),y0(iT-5:iT),'ko-','markersize',10);
                    P5 = plot(nx0(iT-1:iT),ny0(iT-1:iT),'mx-','markersize',10);
                    P6 = plot(tnx(iT-1:iT),tny(iT-1:iT),'bs-','markersize',10);
                else
                    P1 = [];
                    P5 = [];
                    P6 = [];
                end
                P2 = plot(x0(iT),y0(iT),'ro','markersize',10);
                P3 = plot(nx0(iT),ny0(iT),'mx','markersize',10,'linewidth',2);
                P4 = plot(tnx(iT),tny(iT),'bs','markersize',10,'linewidth',2);
                
                
                
                title(sprintf('Frame %d, x(t)=%0.2f (%0.2f)',frame0(iT),x0(iT),time0(iT)),'fontsize',24);
                
                
                if isnan(nx0(iT))
                else
                    pause;
                end
                delete(P1);
                delete(P2);
                delete(P3);
                delete(P4);
                delete(P5);
                delete(P6);
                delete(P7);
            end
        else
            tnx(iT) = nx0(iT);
            tny(iT) = ny0(iT);
            ttx(iT) = tx0(iT);
            tty(iT) = ty0(iT);
        end
    end
end







