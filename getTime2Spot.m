%getTime2Spot
% 2017-07-13 AndyP
% gets time to spot / trail

doTest = false;
threshold = 0.5; % cm (max dist to spot to indicate spot is found)
leavethreshold = 10; % cm (max dist from spot to indicate spot is left)
minT = 0.2; % s (minimum time @ spot to consider spot found)
postSmoothing = 0.1; % 



nSm = ceil(postSmoothing/(1/50));
nM = length(unique(mouse));
nS = max(sess);

initD = nan(nM,nS);
t2s = nan(nM,nS);
minD = nan(nM,nS);
dwellT = nan(nM,nS);
sess0 = nan(nM,nS);
mouse0 = nan(nM,nS);
spotfound = nan(nM,nS);
quadrant = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS & bait==0;
        
        dnT0 = dnT(k0);
        frame0 = frame(k0);
        
        kT = mouseT==iM & sessT==iS;
        
        xT0 = nanmedian(xT1(kT));
        yT0 = nanmedian(yT1(kT));
        
        if xT0 > 1280/2 && yT0 > 1018/2
            quadrant(iM,iS)=1;
        elseif xT0 > 1280/2 && yT0 < 1018/2
            quadrant(iM,iS)=2;
        elseif xT0 < 1280/2 && yT0 > 1018/2
            quadrant(iM,iS)=3;
        elseif xT0 < 1280/2 && yT0 < 1018/2
            quadrant(iM,iS)=4;
        end
    
        %if sum(k0)>0
        if conc(k0)==2
            sess0(iM,iS)=iS;
            mouse0(iM,iS)=iM;
            initD(iM,iS) = dnT0(find(~isnan(dnT0),1,'first'));
            minD(iM,iS) = nanmin(dnT0);
            firstT = find(dnT0<threshold,1,'first');
            if ~isempty(firstT);
                t2s(iM,iS)=frame0(firstT)/50;  % frames -> s
                % find dwell time @ spot
                dnT0temp = dnT0(firstT:end);
                dnT0temp = nanfastsmooth(dnT0temp,nSm,1,2);
                frametemp = frame0(firstT:end);
                lastT = find(dnT0temp>leavethreshold,1,'first');
                if isempty(lastT)
                    lastT = length(dnT0temp);
                end
                
                if doTest
                    F = figure(1); clf; %#ok<UNRCH>
                    plot(frame0,dnT0); hold on;
                    plot(frametemp,dnT0temp,'k');
                    plot(firstT,dnT0(firstT),'r.','markersize',20);
                    plot(lastT+firstT,dnT0temp(lastT),'b.','markersize',20);
                    pause;
                end
                
                dwellT(iM,iS)=(frametemp(lastT)-frame0(firstT))/50;
            else
                t2s(iM,iS)= -20;
            end
            
            spotfound(iM,iS)=dwellT(iM,iS)>minT & t2s(iM,iS)>-1 & ~isnan(t2s(iM,iS));
            
        end
    end
end


disp('Percent of successful trials:');
disp(sum(dwellT(:)>minT & t2s(:)>-1 & ~isnan(t2s(:)))./sum(~isnan(t2s(:))));


