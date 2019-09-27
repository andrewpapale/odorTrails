function [initD,t2s,spotfound,dwellT,quadrant,initO,initOb,spotring,dtraveled]=getTime2Spot2(mouse,sess,dnT,x,y,nx,ny,frame,xT1,yT1)

%getTime2Spot2
% 2017-07-13 AndyP
% 2017-11-28 AndyP same as getTime2SPot but a function
% gets time to spot / trail


threshold = 2; % cm (max dist to spot to indicate spot is found)
leavethreshold = 10; % cm (max dist from spot to indicate spot is left)
minT = 0.5; % s (minimum time @ spot to consider spot found)
postSmoothing = 0.1; % s



nSm = ceil(postSmoothing/(1/50));
nM = max(mouse);
nS = max(sess);

initD = nan(nM,nS);
t2s = nan(nM,nS);
%minD = nan(nM,nS);
dwellT = nan(nM,nS);
%sess0 = nan(nM,nS);
%mouse0 = nan(nM,nS);
%trial0 = nan(nM,nS);
spotfound = nan(nM,nS);
quadrant = nan(nM,nS);
%conc0 = nan(nM,nS);
%bait0 = nan(nM,nS);
% t2s2 = nan(nM,nS);
% t2s3 = nan(nM,nS);
%ALorKP0 = nan(nM,nS);
initO = nan(nM,nS);
initOb = nan(nM,nS);
spotring = nan(nM,nS);
dtraveled = nan(nM,nS);
%nframe = nan(nM,nS);
%meanV = nan(nM,nS);
%meannV = nan(nM,nS);
%xT = nan(nM,nS);
%yT = nan(nM,nS);
iC = 0;
dRing = linspace(0,1000,10);
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS;
        
        
        if sum(k0)>0
            
            iC = iC+1;
            %nframe(iM,iS)=nansum(k0);
            
            dnT0 = dnT(k0);
            frame0 = frame(k0);
            x0 = x(k0);
            y0 = y(k0);
            nx0 = nx(k0);
            ny0 = ny(k0);

            % for shuffled spots
            xT0 = xT1(iC);
            yT0 = yT1(iC);
            
            spotring0 = sqrt((xT0-1280/2).^2+(yT0-1024/2).^2);
            spotring(iM,iS)=find(hist(spotring0,dRing)==1,1,'first');
            
            if xT0<44.8 | xT0>1280-44.8 | yT0<44.8 | yT0>1024-44.8 %#ok<OR2>
                t2s(iM,iS)=-40;
            else
                
                if xT0 > 1280/2 && yT0 > 1024/2
                    quadrant(iM,iS)=2;
                elseif xT0 > 1280/2 && yT0 < 1024/2
                    quadrant(iM,iS)=4;
                elseif xT0 < 1280/2 && yT0 > 1024/2
                    quadrant(iM,iS)=1;
                elseif xT0 < 1280/2 && yT0 < 1024/2
                    quadrant(iM,iS)=3;
                end
                
                %sess0(iM,iS)=iS;
                %mouse0(iM,iS)=iM;
                %trial0(iM,iS)=mode(trial(k0));
                %conc0(iM,iS)=mode(conc(k0));
                %bait0(iM,iS)=mode(bait(k0));
                %ALorKP0(iM,iS)=mode(ALorKP(k0));
                %meanV(iM,iS)=nanmean(V(k0));
                %meannV(iM,iS)=nanmean(nV(k0));
                %xT(iM,iS) = xT0;
                %yT(iM,iS) = yT0;
                %disp(ALorKP0(iM,iS));
                %minD(iM,iS) = nanmin(dnT0);
                firstT = find(dnT0<threshold,1,'first');
                %A = [xT0-x0(1),yT0-y0(1)]./sqrt((yT0-y0(1)).^2+(xT0-x0(1)).^2);
                %B = [nx0(1)-x0(1),ny0(1)-y0(1)]./sqrt((ny0(1)-y0(1)).^2+(nx0(1)-x0(1)).^2);
                %initO(iM,iS)=dot(B,A)*180/pi;
                k2 = find(~isnan(dnT0) & frame0 <= 50,1,'first');
                if ~isempty(k2)
                   A = nanmean(atan2(repmat(yT0,[10,1])-y0(k2:k2+9),repmat(xT0,[10,1])-x0(k2:k2+9)));
                    B = nanmean(atan2(ny0(k2:k2+9)-y0(k2:k2+9),nx0(k2:k2+9)-x0(k2:k2+9)));
                    initO(iM,iS)=angdiff(A,B)*180/pi;
                    initOb(iM,iS)=A*180/pi;
                    initD(iM,iS)=dnT0(k2);
                end
                if ~isempty(firstT)
                    
                    t2s(iM,iS)=frame0(firstT)/50;  % frames -> s
                    
                    if ~isempty(k2)
                        dtraveled(iM,iS) = nansum(sqrt(nx0(k2:firstT).^2+ny0(k2:firstT).^2));
                    end                    
                    
                    % find dwell time @ spot
                    dnT0temp = dnT0(firstT:end);
                    dnT0temp = nanfastsmooth(dnT0temp,nSm,1,2);
                    frametemp = frame0(firstT:end);
                    lastT = find(dnT0temp>leavethreshold,1,'first');
                    if isempty(lastT)
                        lastT = length(dnT0temp);
                    end
                    
%                     % find 2nd time to spot (if any)
%                     secondT = find(dnT0<threshold & frame0 > lastT+firstT,1,'first');
%                     if ~isempty(secondT)
%                         t2s2(iM,iS)=frame0(secondT)/50;
%                     else
%                         t2s2(iM,iS) = -20;
%                     end
                    dwellT(iM,iS)=(frametemp(lastT)-frame0(firstT))/50;
                else
                    t2s(iM,iS)= -20;
                end
                
                spotfound(iM,iS)=dwellT(iM,iS)>minT & t2s(iM,iS)>-1 & ~isnan(t2s(iM,iS));
                
            end
        end
    end
    disp(iM);
end


disp('Percent of successful trials:');
disp(sum(dwellT(:)>minT & t2s(:)>-1 & ~isnan(t2s(:)))./sum(~isnan(t2s(:)) & t2s(:)~=-40));


