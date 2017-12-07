% get orientation from spot

doTest = false;

nM = max(mouse);
nS = max(sess);

orient = nan(size(x));
borient = nan(size(x));
bsorient = nan(size(x));
%dO = nan(size(x));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            kT = mouseT==iM & sessT==iS;
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
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
            
            A = atan2(yT-y0,xT-x0);
            B = atan2(ny0-y0,nx0-x0);
            %orient(k) = A;
            orient(k) = angdiff(A,B)*180/pi;
%             orient1 = nan(size(x0));
%             k1 = ~isnan(A) & ~isnan(B);
%             orient0=acos(dot(A(k1),B(k1),2))*180/pi;
%             orient1(k1)=orient0;
%             orient(k)=orient1;
            
            [ th_head,~,borient0] = nose_deflect(x0,y0,nx0,ny0);
            borient(k)=borient0*180/pi;
            bsorient(k)=angdiff(th_head,borient(k)*pi/180)*180/pi;
            
            %dO(k) = foaw_diff(orient(k),1/50,50,0.5,0.1);
            
            
            if doTest
                F = figure(1); clf; 
                
                subplot(1,2,1);
                plot(x0,y0,'k.'); hold on;
                plot(nx0,ny0,'b.');
                for iP=1:length(x0)
                    line([x0(iP) nx0(iP)],[y0(iP) ny0(iP)],'color','k');
                end
                xA = cos(A);
                yA = sin(A);
                quiver(x0,y0,xA,yA);
                xB = cos(B);
                yB = sin(B);
                quiver(x0,y0,xB,yB,'color','c');
                plot(xT,yT,'r.','markersize',30);
                
                subplot(1,2,2);
                plot(x0,y0,'k.'); hold on;
                plot(nx0,ny0,'b.');
                for iP=1:length(x0)
                    line([x0(iP) nx0(iP)],[y0(iP) ny0(iP)],'color','k');
                end
                xC = cos((orient(k)*pi/180));
                yC = sin((orient(k)*pi/180));
                quiver(x0,y0,xC,yC);
                plot(xT,yT,'r.','markersize',30);
                pause;
            end
            
        end
    end
    disp(iM);
end