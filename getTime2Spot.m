%getTime2Spot
% 2017-07-13 AndyP
% gets time to spot / trail


doTest = false;
threshold = 2; % cm (max dist to spot to indicate spot is found)
leavethreshold = 10; % cm (max dist from spot to indicate spot is left)
minT = 0; % s (minimum time @ spot to consider spot found)
postSmoothing = 0.1; % s



nSm = ceil(postSmoothing/(1/50));
nM = length(unique(mouse));
nS = max(sess);

initD = nan(nM,nS);
initDb = nan(nM,nS);
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
ALorKP0 = nan(nM,nS);
initO = nan(nM,nS);
nframe = nan(nM,nS);
meanV = nan(nM,nS);
meannV = nan(nM,nS);
medianV = nan(nM,nS);
mediannV = nan(nM,nS);
meandphi = nan(nM,nS);
meandphiV = nan(nM,nS);
fractoward = nan(nM,nS);
fracedge = nan(nM,nS);
ddV = nan(nM,nS);
ddnV = nan(nM,nS);
xT = nan(nM,nS);
yT = nan(nM,nS);
spotring = nan(nM,nS);

initX = nan(nM,nS);
initY = nan(nM,nS);
initnX = nan(nM,nS);
initnY = nan(nM,nS);
Squad = nan(nM,nS);

texplo = nan(nM,nS);


V1 = V;
nV1 = nV;

dRing = linspace(0,1000,10);

for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS;
        kedge = mouse==iM & sess==iS & edge;
        
        if sum(k0)>0
            
            V1(kedge)=nan;
            nV1(kedge)=nan;
            nframe(iM,iS)=nansum(k0);
            
            dnT0 = dnT(k0);
            frame0 = frame(k0);
            x0 = x(k0);
            y0 = y(k0);
            nx0 = nx(k0);
            ny0 = ny(k0);
            
            kT = mouseT==iM & sessT==iS;
            
            xT0 = nanmedian(xT1(kT));
            yT0 = nanmedian(yT1(kT));
            
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
                
                sess0(iM,iS)=iS;
                mouse0(iM,iS)=iM;
                trial0(iM,iS)=mode(trial(k0));
                conc0(iM,iS)=mode(conc(k0));
                bait0(iM,iS)=mode(bait(k0));
                ALorKP0(iM,iS)=mode(ALorKP(k0));
%                 meanV(iM,iS)=nanmean(V(k0));
%                 medianV(iM,iS)=nanmedian(V(k0));
%                 meannV(iM,iS)=nanmean(nV(k0));
%                 mediannV(iM,iS)=nanmedian(nV(k0));
%                 k1 = k0 & nV>0.01;
%                 meandphi(iM,iS)=nanmean(abs(nidphi(k1))./nV(k1));
                xT(iM,iS) = xT0;
                yT(iM,iS) = yT0;
                %disp(ALorKP0(iM,iS));
                initD(iM,iS) = dnT0(find(~isnan(dnT0),1,'first'));
                initDb(iM,iS) = dnT(find(~isnan(dnT),1,'first'));
                
                
                minD(iM,iS) = nanmin(dnT0);
                firstT = find(dnT0<threshold,1,'first');
                %A = [xT0-x0(1),yT0-y0(1)]./sqrt((yT0-y0(1)).^2+(xT0-x0(1)).^2);
                %B = [nx0(1)-x0(1),ny0(1)-y0(1)]./sqrt((ny0(1)-y0(1)).^2+(nx0(1)-x0(1)).^2);
                %initO(iM,iS)=dot(B,A)*180/pi;
                
                xI = x0(V1(k0) > 0.5);
                yI = y0(V1(k0) > 0.5);
                nxI = nx0(V1(k0)> 0.5);
                nyI = ny0(V1(k0)> 0.5);
                xI = xI(1);
                yI = yI(1);
                nxI = nxI(1);
                nyI = nyI(1);
                
                initX(iM,iS) = xI;
                initY(iM,iS) = yI;
                initnX(iM,iS)=nxI;
                initnY(iM,iS)=nyI;
                
                if initX(iM,iS) > 1280/2 && initY(iM,iS) > 1024/2
                    Squad(iM,iS)=2;
                elseif initX(iM,iS) > 1280/2 && initY(iM,iS) < 1024/2
                    Squad(iM,iS)=4;
                elseif initX(iM,iS) < 1280/2 && initY(iM,iS) > 1024/2
                    Squad(iM,iS)=1;
                elseif initX(iM,iS) < 1280/2 && initY(iM,iS) < 1024/2
                    Squad(iM,iS)=3;
                end
                
                A = atan2(yT0-yI,xT0-xI);
                B = atan2(nyI-yI,nxI-xI);
                initO(iM,iS)= B*180/pi;
                if ~isempty(firstT)
                    
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
                    if ~isempty(secondT)
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
                        if ~isempty(thirdT)
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
                        plot(frame0/50,dnT0,'linewidth',2); hold on;
                        plot(frametemp/50,dnT0temp,'k','linewidth',2);
                        plot(firstT/50,dnT0(firstT),'r.','markersize',50);
                        plot((lastT+firstT)/50,dnT0temp(lastT),'b.','markersize',50);
                        if ~isnan(secondT)
                            plot(secondT/50,dnT0(secondT),'r.','markersize',50);
                            plot((lastT2+secondT)/50,dnT0temp2(lastT2),'b.','markersize',50);
                            if ~isnan(thirdT)
                                plot(thirdT/50,dnT0(thirdT),'r.','markersize',50);
                                plot((lastT3+thirdT)/50,dnT0temp3(lastT3),'b.','markersize',50);
                            end
                        end
                        pause;
                    end
                    
                    dwellT(iM,iS)=(frametemp(lastT)-frame0(firstT))/50;
                else
                    t2s(iM,iS)= -20;
                end
                
                tempvec = zeros(sum(k0),1);
                tempedge = find(k0==1 & kedge==1);
                firstframe = find(k0==1,1,'first');
                tempedge = tempedge-firstframe+1;
                if t2s(iM,iS)>0
                    tempvec(1:nanmin(round(t2s(iM,iS)*50)+25,length(tempvec)))=1;
                else
                    dummytime = nanmean(t2s(iM,:)*50)+25;
                    tempvec(1:nanmin(round(dummytime),length(tempvec)))=1;
                end
                texplo(iM,iS) = sum(tempvec)/50;
                tempedge(tempedge>length(tempvec))=[];
                tempvec(tempedge)=0;
                V0 = V1(k0);
                nV0 = nV1(k0);
                nidphi0 = nidphi(k0);
                meanV(iM,iS)=nanmean(V0(tempvec==1));
                medianV(iM,iS)=nanmedian(V0(tempvec==1));
                meannV(iM,iS)=nanmean(nV0(tempvec==1));
                mediannV(iM,iS)=nanmedian(nV0(tempvec==1));
                k1 = tempvec & nV0>0.01;
                meandphi(iM,iS)=nanmean(abs(nidphi0(k1))./nV0(k1));
                meandphiV(iM,iS) = nanmean(abs(nidphi0(k1)./nV0(k1))./V0(k1));
                fractoward(iM,iS) = nansum(dnT(tempvec==1)>dT(tempvec==1))./sum(tempvec==1);
                fracedge(iM,iS) = length(tempedge)./length(tempvec);
                ddnV(iM,iS) = nanmean(foaw_diff(nV0(tempvec==1),1/50,50,0.5,0.2));
                ddV(iM,iS) = nanmean(foaw_diff(V0(tempvec==1),1/50,50,0.5,0.2));
                
                spotfound(iM,iS)=dwellT(iM,iS)>minT & t2s(iM,iS)>-1 & ~isnan(t2s(iM,iS));
                
            end
        end
    end
    disp(iM);
end


disp('Percent of successful trials:');
disp(sum(dwellT(:)>minT & t2s(:)>-1 & ~isnan(t2s(:)))./sum(~isnan(t2s(:)) & t2s(:)~=-40));

