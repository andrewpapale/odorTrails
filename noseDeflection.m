% test nose deflection from centerline
% 2018-08-11 AndyP, fixed 2 errors, line 53, line 81, distance computation,
% important!
%
doTest = false;

nM = max(mouse);
nS = max(sess);
Dline = nan(size(mouse));
% Dline2 = nan(size(mouse));
% whichSide = nan(size(mouse));
% whichSideSpot = nan(size(mouse));
dwSide = nan(size(mouse));
slope = nan(size(mouse));
theta = nan(size(mouse));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        %kT = mouseT==iM & sessT==iS;
        if sum(k)>0
            
%             cx0 = nanmedian(xT1(kT));
%             cy0 = nanmedian(yT1(kT));
            
            cx0 = cx(iS);
            cy0 = cy(iS);

            %edge0 = edge(k);
            %isstopped0 = isstopped(k);
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            %nx0(edge0 & isstopped0)=nan;
            %ny0(edge0 & isstopped0)=nan;
            
            if doTest
                figure(1); clf;
                plot(x0,y0,'.','markersize',1,'color',[0.5 0.5 0.5]);
                hold on;
                viscircles([cx0 cy0],2);
            end
            
            dx = foaw_diff(x0,1/50,50,0.2,0.1);
            dy = foaw_diff(y0,1/50,50,0.2,0.1);
            dy0 = foaw_diff(y0,1/50,50,0.2,0.1);
            dx0 = foaw_diff(x0,1/50,50,0.2,0.1);
            m0 = dy0./dx0;
            
            minX = 1;
            maxX = 1300;
            
            
            nT = length(x0);
            Dline0 = nan(nT,1);
            theta0 = nan(nT,1);
            
            for iT=1:nT
                result = computeline([x0(iT), y0(iT)],m0(iT), [minX maxX]);
                nP = length(result);
                x1 = result{1}(:,1);
                y1 = result{1}(:,2);
                x2 = result{nP}(:,1);
                y2 = result{nP}(:,2);
                
                Dline0(iT) = ((y2-y1)*nx0(iT)-(x2-x1)*ny0(iT)+x2*y1-y2*x1)./sqrt((y2-y1).^2+(x2-x1).^2)./11.2;
                dxline = sqrt((nx0(iT)-x0(iT)).^2+(ny0(iT)-x0(iT)).^2)./11.2;
                theta0(iT) = acosd(dxline/sqrt(dxline.^2+Dline0(iT).^2));
                
                if doTest
                    %P7 = plot(x1,y1,'r.-');
                    P7 = line([x1 x2],[y1 y2],'color','r');
                    P8 = plot(x3(ix),y3(ix),'g.','markersize',10);
                    if iT>5
                        P1 = plot(x0((iT)-3:(iT)),y0((iT)-3:iT),'ko-','markersize',10);
                        P5 = plot(nx0((iT)-1:(iT)),ny0((iT)-1:iT),'mx-','markersize',10);
                        P6 = [];
                    else
                        P1 = [];
                        P5 = [];
                        P6 = [];
                    end
                    P2 = plot(x0(iT),y0(iT),'ro','markersize',10);
                    
                    %                         switch whichSide0(totest(iT))
                    %                             case 0
                    P3 = plot(nx0(iT),ny0(iT),'rx','markersize',10,'linewidth',2);
                    %                             case 1
                    %                               P3 = plot(nx0(totest(iT)),ny0(totest(iT)),'mx','markersize',10,'linewidth',2);
                    %                         end
                    P4 = [];
                    
                    set(gca,'XLim',[1 1300],'Ylim',[1 1050]);
                    
                    title(sprintf('%0.1d',Dline0),'fontsize',24);
                    
                    
                    if isnan(Dline0)
                    else
                        %disp(Dline0(totest(iT)));
                        %disp(Dline3(totest(iT)));
                        %disp(Dline0(totest(iT))-Dline3(totest(iT)));
                        pause;
                    end
                    delete(P1);
                    delete(P2);
                    delete(P3);
                    delete(P4);
                    delete(P5);
                    delete(P6);
                    delete(P7);
                    delete(P8);
                end
               
            end
            Dline(k)=Dline0;
            theta(k)=theta0;
        end
        disp(iS);
    end
    disp(iM);
end