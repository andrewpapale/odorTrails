%getTime2Spot3
% 2017-07-13 AndyP
% gets time to spot / trail
thrVec = 3;
%thrVec = [0.1 0.25 0.33 0.5 0.66 0.75 1 1.25 1.5 2 3 4 5];
%thrVec = 0.6:0.02:0.8; % sensitivity analysis
%thrVec = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 1.75 2 3 4 5];
nR = length(thrVec);
initD = nan(nM,nS,nR);
initDb = nan(nM,nS,nR);
t2s = nan(nM,nS,nR);
minD = nan(nM,nS,nR);
dwellT = nan(nM,nS,nR);
sess0 = nan(nM,nS,nR);
mouse0 = nan(nM,nS,nR);
trial0 = nan(nM,nS,nR);
spotfound = nan(nM,nS,nR);
quadrant = nan(nM,nS,nR);
conc0 = nan(nM,nS,nR);
bait0 = nan(nM,nS,nR);
t2s2 = nan(nM,nS,nR);
t2s3 = nan(nM,nS,nR);
ALorKP0 = nan(nM,nS,nR);
initO = nan(nM,nS,nR);
initObn = nan(nM,nS,nR);
nframe = nan(nM,nS,nR);
meanV = nan(nM,nS,nR);
meannV = nan(nM,nS,nR);
medianV = nan(nM,nS,nR);
mediannV = nan(nM,nS,nR);
meandphi = nan(nM,nS,nR);
idphi = nan(nM,nS,nR);

meandphiV = nan(nM,nS,nR);
fractoward = nan(nM,nS,nR);
fracedge = nan(nM,nS,nR);
ddV = nan(nM,nS,nR);
ddnV = nan(nM,nS,nR);
xT = nan(nM,nS,nR);
yT = nan(nM,nS,nR);
spotring = nan(nM,nS,nR);
dtraveled = nan(nM,nS,nR);
spotfound2 = nan(nM,nS,nR);
angcast0 = nan(nM,nS,nR);
tstopped = nan(nM,nS,nR);

initX = nan(nM,nS,nR);
initY = nan(nM,nS,nR);
initnX = nan(nM,nS,nR);
initnY = nan(nM,nS,nR);
Squad = nan(nM,nS,nR);
ThrV = nan(nM,nS,nR);
V1 = V;
nV1 = nV;

texplo = nan(nM,nS,nR);

weight0 = nan(nM,nS,nR);
pweight0 = nan(nM,nS,nR);

for iRc=1:nR
doTest = false;
threshold = thrVec(iRc); %0.72; %thrVec(iR); % cm (max dist to spot to indicate spot is found)
leavethreshold = 10; % cm (max dist from spot to indicate spot is left)
minT = 0; %thrVec(iR); % s (minimum time @ spot to consider spot found)
postSmoothing = 0.1; % s



nSm = ceil(postSmoothing/(1/50));
nM = length(unique(mouse));
nS = max(sess);



dRing = linspace(0,1000,5);
lastt2s = 1500;
iC = 1;
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS & V< 100 & nV < 100;
        kedge = k0 & edge;
        kstop0 = k0 & isstopped & V < 100;
        
        if sum(k0)>0
            
            ThrV(iM,iS,iRc)=iRc;
            V1(kedge)=nan;
            nV1(kedge)=nan;
            nframe(iM,iS,iRc)=nansum(k0);
            
            dnT0 = dnT(k0);
            dT0 = dT(k0);
            frame0 = frame(k0);
            x0 = x(k0);
            y0 = y(k0);
            nx0 = nx(k0);
            ny0 = ny(k0);
            orient0 = orient(k0);
            
            kT = mouseT==iM & sessT==iS;
            
            xT0 = cx0(iM,iS);
            yT0 = cy0(iM,iS);
            
%             spotring0 = sqrt((xT0-1280/2).^2+(yT0-1024/2).^2);
%             spotring(iM,iS,iRc)=find(hist(spotring0,dRing)==1,1,'first');
%             
            %weight0(iM,iS,iRc) = weight(iC);
            %pweight0(iM,iS,iRc) = perweight(iC);
            iC = iC+1;
            
            if xT0<44.8 | xT0>1280-44.8 | yT0<44.8 | yT0>1024-44.8 %#ok<OR2>
                t2s(iM,iS,iRc)=-40;
            else
                
                if xT0 > 1280/2 && yT0 > 1024/2
                    quadrant(iM,iS,iRc)=2;
                elseif xT0 > 1280/2 && yT0 < 1024/2
                    quadrant(iM,iS,iRc)=4;
                elseif xT0 < 1280/2 && yT0 > 1024/2
                    quadrant(iM,iS,iRc)=1;
                elseif xT0 < 1280/2 && yT0 < 1024/2
                    quadrant(iM,iS,iRc)=3;
                end
                
                sess0(iM,iS,iRc)=iS;
                mouse0(iM,iS,iRc)=iM;
%                 if ~isempty(trial)
%                     trial0(iM,iS,iRc)=mode(trial(k0));
%                 end
%                 if ~isempty(conc)
%                     conc0(iM,iS,iRc)=mode(conc(k0));
%                 end
%                 if ~isempty(bait)
%                     bait0(iM,iS,iRc)=mode(bait(k0));
%                 end
%                 if ~isempty(ALorKP)
%                     ALorKP0(iM,iS,iRc)=mode(ALorKP(k0));
%                 end
                
%                 meanV(iM,iS,iRc)=nanmean(V(k0));
%                 medianV(iM,iS,iRc)=nanmedian(V(k0));
%                 meannV(iM,iS,iRc)=nanmean(nV(k0));
%                 mediannV(iM,iS,iRc)=nanmedian(nV(k0));
%                 k1 = k0 & nV>0.01;
%                 meandphi(iM,iS,iRc)=nanmean(abs(nidphi(k1))./nV(k1));
                xT(iM,iS,iRc) = xT0;
                yT(iM,iS,iRc) = yT0;
                %disp(ALorKP0(iM,iS,iRc));
                %k2 = find(~isnan(dnT0) & frame0 <= 50,1,'first');
                initD(iM,iS,iRc) = sqrt((cx0(iM,iS)-575).^2+(cy0(iM,iS)-15).^2);
                initDb(iM,iS,iRc) = dT0(1);
                initX(iM,iS,iRc) = x0(1);
                initY(iM,iS,iRc) = y0(1);
                initnX(iM,iS,iRc) = nx0(1);
                initnY(iM,iS,iRc) = ny0(1);
                initO(iM,iS,iRc)= orient0(1);
                %initObn(iM,iS,iRc)=B*180/pi;
               if initX(iM,iS,iRc) > 1280/2 && initY(iM,iS,iRc) > 1024/2
                    Squad(iM,iS,iRc)=2;
                elseif initX(iM,iS,iRc) > 1280/2 && initY(iM,iS,iRc) < 1024/2
                    Squad(iM,iS,iRc)=4;
                elseif initX(iM,iS,iRc) < 1280/2 && initY(iM,iS,iRc) > 1024/2
                    Squad(iM,iS,iRc)=1;
                elseif initX(iM,iS,iRc) < 1280/2 && initY(iM,iS,iRc) < 1024/2
                    Squad(iM,iS,iRc)=3;
                end
                
                
                minD(iM,iS,iRc) = nanmin(dnT(k0));
                firstT = find(dnT0<threshold,1,'first');
                %A = [xT0-x0(1),yT0-y0(1)]./sqrt((yT0-y0(1)).^2+(xT0-x0(1)).^2);
                %B = [nx0(1)-x0(1),ny0(1)-y0(1)]./sqrt((ny0(1)-y0(1)).^2+(nx0(1)-x0(1)).^2);
                %initO(iM,iS,iRc)=dot(B,A)*180/pi;
                if ~isempty(firstT)
                    
                    if ~isempty(k2)
                        k3 = k2:firstT;
                        k3(isnan(k3))=[];
                        dtraveled(iM,iS,iRc) = nansum(sqrt(diff(nx0(k3)).^2+diff(ny0(k3)).^2))./11.2;
                    end 
                    
                    t2s(iM,iS,iRc)=frame0(firstT)/50;  % frames -> s
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
                        t2s2(iM,iS,iRc)=frame0(secondT)/50;
                        dnT0temp2 = dnT0(secondT:end);
                        %dnT0temp2 = nanfastsmooth(dnT0temp2,nSm,1,2);
                        frametemp2 = frame0(secondT:end);
                        lastT2 = find(dnT0temp2>leavethreshold,1,'first');
                        if isempty(lastT2)
                            lastT2 = length(dnT0temp2);
                        end
                       
                        % find 3rd time to spot (if any)
                        thirdT = find(dnT0<threshold  & frame0 > lastT2+secondT,1,'first');
                        if ~isempty(thirdT)
                            t2s3(iM,iS,iRc)=frame0(thirdT)/50;
                            dnT0temp3 = dnT0(thirdT:end);
                            %dnT0temp3 = nanfastsmooth(dnT0temp3,nSm,1,2);
                            frametemp3 = frame0(thirdT:end);
                            lastT3 = find(dnT0temp3>leavethreshold,1,'first');
                            if isempty(lastT3)
                                lastT3 = length(dnT0temp3);
                            end
                        else
                            t2s3(iM,iS,iRc) = -20;
                        end
                    else
                        secondT = nan;
                        t2s2(iM,iS,iRc) = -20;
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
                    
                    dwellT(iM,iS,iRc)=(frametemp(lastT)-frame0(firstT))/50;
                    
%                     if ~isempty(lastT2) && ~isnan(secondT)
%                         dwellT2 =(frametemp2(lastT2)-frame0(secondT))/50;
%                     else
%                         dwellT2 = nan;
%                     end
                else
                    t2s(iM,iS,iRc)= -20;
                end
                
                tempvec = zeros(sum(k0),1);
                tempedge = find(k0==1 & kedge==1);
                tempstop = find(k0==1 & kstop0==1);
                firstframe = find(k0==1,1,'first');
                tempedge = tempedge-firstframe+1;
                tempstop = tempstop-firstframe+1;
                if t2s(iM,iS,iRc)>0
                    tempvec(1:nanmin(round(t2s(iM,iS,iRc)*50),length(tempvec)))=1;
                    if t2s(iM,iS,iRc)*50+25 > length(tempvec)
                        endflag = 1;
                    else
                        endflag = 0;
                    end
                    lastt2s = round(t2s(iM,iS,iRc));
                else
                    dummytime = lastt2s*50;
                    
                    tempvec(1:nanmin(round(dummytime),length(tempvec)))=1;
                    if round(dummytime) >length(tempvec)
                        endflag = 1;
                    else
                        endflag = 0;
                    end
                end
%                 if any(tempedge>length(tempvec))
%                     keyboard;
%                 end
                tempedge(tempedge>length(tempvec))=[]; % ok, none here
                tempstop(tempstop>length(tempvec))=[];
                tempvec(tempedge)=0;
                if ~endflag
                    texplo(iM,iS,iRc)=(sum(tempvec==1))/50;
                else
                  texplo(iM,iS,iRc)=sum(tempvec==1)/50;
                end
                tempvec(tempstop)=0;
                V0 = V1(k0);
                nV0 = nV1(k0);
                nidphi0 = ndphi(k0);
                if sum(tempvec)>25
                    meanV(iM,iS,iRc)=nanmean(V0(tempvec==1));
                    medianV(iM,iS,iRc)=nanmedian(V0(tempvec==1));
                    meannV(iM,iS,iRc)=nanmean(nV0(tempvec==1));
                    mediannV(iM,iS,iRc)=nanmedian(nV0(tempvec==1));
                    k1 = tempvec & nV0>0.01 & V0 > 0.01;
                    log10dphi0 = logVar(abs(nidphi0)./nV0);
                    meandphi(iM,iS,iRc)=nanmean(log10dphi0(k1));
                    meandphiV(iM,iS,iRc) = nanmean(log10dphi0(k1)./V0(k1));
                    idphi(iM,iS,iRc) = nansum(abs(nidphi0(k1))./nV0(k1));
                    fractoward(iM,iS,iRc) = nansum(dnT(tempvec==1)>dT(tempvec==1))./sum(tempvec==1);
                    fracedge(iM,iS,iRc) = length(tempedge)./length(tempvec);
                    ddnV(iM,iS,iRc) = nanmean(foaw_diff(nV0(tempvec==1),1/50,50,0.2,0.1));
                    ddV(iM,iS,iRc) = nanmean(foaw_diff(V0(tempvec==1),1/50,50,0.2,0.1));
                    %angcast0(iM,iS,iRc)=nansum(angcast(tempvec==1))./sum(tempvec==1);
                    tstopped(iM,iS,iRc)=length(tempstop)/50;
                end
                spotfound(iM,iS,iRc)=dwellT(iM,iS,iRc)>minT & t2s(iM,iS,iRc)>0;
                %spotfound2(iM,iS,iRc)=dwellT2>minT & t2s2(iM,iS,iRc)>0;
            end
        end
    end
    disp(iM);
end


%disp('Percent of successful trials:');
%disp(sum(dwellT(:)>minT & t2s(:)>-1 & ~isnan(t2s(:)))./sum(~isnan(t2s(:)) & t2s(:)~=-40));

iRc = iRc+1;
disp(iRc);

end
