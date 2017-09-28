%getTime2Spot
% 2017-07-13 AndyP
% gets time to spot / trail

doTest = false;
threshold = 1; % cm (max dist to spot to indicate spot is found)
leavethreshold = 10; % cm (max dist from spot to indicate spot is left)
minT = 0.2; % s (minimum time @ spot to consider spot found)
postSmoothing = 0.1; % s



nSm = ceil(postSmoothing/(1/50));
nM = length(unique(mouse));
nS = max(sess);

initD = nan(nM,nS);
t2s = nan(nM,nS);
minD = nan(nM,nS);
dwellT = nan(nM,nS);
sess0 = nan(nM,nS);
mouse0 = nan(nM,nS);
trial0 = nan(nM,nS);
spotfound = nan(nM,nS);
quadrant = nan(nM,nS);
conc0 = nan(nM,nS);
bait0 = nan(nM,nS);
t2s2 = nan(nM,nS);
t2s3 = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS;
        
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
        
        k1 = k0 & frame==1;
        
        if sum(k0)>0
            sess0(iM,iS)=iS;
            mouse0(iM,iS)=iM;
            trial0(iM,iS)=trial(k1);
            conc0(iM,iS)=conc(k1);
            bait0(iM,iS)=bait(k1);
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
                
                % find 2nd time to spot (if any)
                secondT = find(dnT0<threshold & frame0 > lastT+firstT,1,'first');
                if ~isempty(secondT);
                    t2s2(iM,iS)=frame0(secondT)/50;
                    dnT0temp2 = dnT0(secondT:end);
                    dnT0temp2 = nanfastsmooth(dnT0temp2,nSm,1,2);
                    frametemp2 = frame0(secondT:end);
                    lastT2 = find(dnT0temp2>leavethreshold,1,'first');
                    if isempty(lastT2)
                        lastT2 = length(dnT0temp2);
                    end
                    
                    % find 3rd time to spot (if any)
                    thirdT = find(dnT0<threshold  & frame0 > lastT2+secondT,1,'first');
                    if ~isempty(thirdT);
                        t2s3(iM,iS)=frame0(thirdT)/50;
                        dnT0temp3 = dnT0(thirdT:end);
                        dnT0temp3 = nanfastsmooth(dnT0temp3,nSm,1,2);
                        frametemp3 = frame0(thirdT:end);
                        lastT3 = find(dnT0temp3>leavethreshold,1,'first');
                        if isempty(lastT3)
                            lastT3 = length(dnT0temp3);
                        end
                    else
                        t2s3(iM,iS) = -20;
                    end
                else
                    secondT = nan;
                    t2s2(iM,iS) = -20;
                end
                
                if doTest
                    F = figure(1); clf; %#ok<UNRCH>
                    plot(frame0,dnT0); hold on;
                    plot(frametemp,dnT0temp,'k');
                    plot(firstT,dnT0(firstT),'r.','markersize',20);
                    plot(lastT+firstT,dnT0temp(lastT),'b.','markersize',20);
                    if ~isnan(secondT)
                        plot(secondT,dnT0(secondT),'r.','markersize',20);
                        plot(lastT2+secondT,dnT0temp2(lastT2),'b.','markersize',20);
                        if ~isnan(thirdT)
                            plot(thirdT,dnT0(thirdT),'r.','markersize',20);
                            plot(lastT3+thirdT,dnT0temp3(lastT3),'b.','markersize',20);
                        end
                    end
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


