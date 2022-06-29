% get orientation from spot

doTest = false;

nM = max(mouse);
nS = max(sess);

orient = [];
% borient = nan(size(x));
% norient = nan(size(x));
torient = [];
forient = [];
%dO = nan(size(x));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            
            xT = cx0(iM,iS);
            yT = cy0(iM,iS);
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            %             nearX0 = nearX(k);
            %             nearY0 = nearY(k);
            %
            %             %// Index array for factor
            %             x1 = 1:numel(x0);
            %             %// Indices of NaNs
            %             t2 = find(~isnan(x0));
            %
            %             %// Replace NaNs with the closest non-NaNs
            %             x0 = interp1(x1(t2),x0(t2),x1,'linear')';
            %             y0 = interp1(x1(t2),y0(t2),x1,'linear')';
            %
            %             %// Index array for factor
            %             nx1 = 1:numel(nx0);
            %             %// Indices of NaNs
            %             t2 = find(~isnan(nx0));
            %
            %             %// Replace NaNs with the closest non-NaNs
            %             nx0 = interp1(nx1(t2),nx0(t2),nx1,'linear')';
            %             ny0 = interp1(nx1(t2),ny0(t2),nx1,'linear')';
            
            BSx = xT-x0;
            BSy = yT-y0;
            BNx = nx0-x0;
            BNy = ny0-y0;
            BFx = 575-x0;
            BFy = 15-y0;
            
            %orient0 = wrapTo180(atan2d(BNy,BNx)-atan2d(BSy,BSx));
            %orient = cat(1,orient,orient0);
            
            %forient0 = wrapTo180(atan2d(BNy,BNx)-atan2d(BFy,BFx));
            %forient = cat(1,forient,forient0);
            
            
%             A = atan2(yT-y0,xT-x0);
%             B = atan2(ny0-y0,nx0-x0);
%             orient(k)=angdiff(A,B)*180/pi;
            torient0 = atan2d(BNy,BNx);
            torient = cat(1,torient,torient0);
            
%             Bx = cat(1,nan,diff(x0));
%             By = cat(1,nan,diff(y0));
%             Nx = cat(1,nan,diff(nx0));
%             Ny = cat(1,nan,diff(ny0));
%             NSx = xT-nx0;
%             NSy = yT-ny0;
            
%             k0 = ~isnan(Bx) & ~isnan(BSx);
%             Ik = find(k0==1);
%             firstk = find(k==1,1,'first');
%             borient0 = nan(sum(k),1);
%             A = sqrt(BSx.^2+BSy.^2);
%             B = sqrt(Bx.^2+By.^2);
%             borient0(k0)= acosd(dot([BSx(k0),BSy(k0)]',[Bx(k0),By(k0)]')'./(A(k0).*B(k0)));
%             borient(Ik+firstk)=borient0(k0);
%  
%             k0 = ~isnan(Nx) & ~isnan(NSx);
%             Ik = find(k0==1);
%             firstk = find(k==1,1,'first');
%             norient0 = nan(sum(k),1);
%             A = sqrt(NSx.^2+NSy.^2);
%             B = sqrt(Nx.^2+Ny.^2);
%             norient0(k0)= acosd(dot([NSx(k0),NSy(k0)]',[Nx(k0),Ny(k0)]')'./(A(k0).*B(k0)));
%             norient(Ik+firstk)=norient0(k0);
            
            
            if doTest
                F = figure(1); clf;
                
                nT = length(x0);
                xC = cos((orient(k)*pi/180));
                yC = sin((orient(k)*pi/180));
                xD = cos((borient(k)*pi/180));
                yD = sin((borient(k)*pi/180));
                for iT=1:nT
                    plot(x0(iT),y0(iT),'m.'); hold on;
                    plot(nx0(iT),ny0(iT),'k.');
                    xA = cos(A);
                    yA = sin(A);
                    xB = cos(B);
                    yB = sin(B);
                    Qa = quiver(x0(iT),y0(iT),xA(iT)*100,yA(iT)*100,'color','r');
                    Qb = quiver(x0(iT),y0(iT),xB(iT)*100,yB(iT)*100,'color','c');
                    Qi = quiver(x0(iT),y0(iT),xC(iT)*100,yC(iT)*100,'color','g');
                    Qd = quiver(x0(iT),y0(iT),xD(iT)*100,yD(iT)*100,'color','b','linewidth',2);
                    plot(xT,yT,'r.','markersize',20);
                    set(gca,'YLim',[0 1050],'XLim',[0 1250]);
                    pause(0.1);
                    delete(Qi); delete(Qa); delete(Qb); delete(Qd);
                end
                %                 subplot(1,2,1);
                %                 scatter(x0,y0,20,1:length(x0),'filled'); hold on;
                %                 scatter(nx0,ny0,20,1:length(nx0),'filled');
                %                 for iP=1:length(x0)
                %                     line([x0(iP) nx0(iP)],[y0(iP) ny0(iP)],'color','k');
                %                 end
                %                 xA = cos(A);
                %                 yA = sin(A);
                %                 quiver(x0,y0,xA,yA);
                %                 xB = cos(B);
                %                 yB = sin(B);
                %                 quiver(x0,y0,xB,yB,'color','c');
                %                 plot(xT,yT,'r.','markersize',30);
                %                 title('Vectors from body to spot and body to nose');
                %
                %                 subplot(1,2,2);
                %                 scatter(x0,y0,20,1:length(x0),'filled'); hold on;
                %                 scatter(nx0,ny0,20,1:length(nx0),'filled');
                %                 for iP=1:length(x0)
                %                     line([x0(iP) nx0(iP)],[y0(iP) ny0(iP)],'color','k');
                %                 end
                %                 xC = cos((orient(k)*pi/180));
                %                 yC = sin((orient(k)*pi/180));
                %                 quiver(x0,y0,xC,yC);
                %                 plot(xT,yT,'r.','markersize',30);
                %                 title('Vector to Spot');
                %                 pause;
            end
            
        end
    end
    disp(iM);
end