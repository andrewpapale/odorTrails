mouse2 = [];
welchSpot = [];
welchCont = [];
welchFar = [];
kSpotS = [];
kSpotF = [];
kStartC = [];
kStopC = [];
dntspot = [];
dntcont = [];
freqs = linspace(0.1,10,100);
spotfound = [];

segLen = 500;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        k2 = nan;
        if sum(k)>0
            nearX0 = nearX(k);
            nearY0 = nearY(k);
            edge0 = ~(nearX0>44.8 & nearY0>44.8 & nearX0<1280-44.8 & nearY0<1024-44.8);
            
            welch0 = nan(size(freqs));
            welch1 = nan(size(freqs));
            welch2 = nan(size(freqs));
            
            dnT0 = dnT(k);
            orient0 = Iorient(k);
            orient0(edge0)=nan;
            if t2s(iM,iS)>0 % found spot
                k2 = round(t2s(iM,iS)*50);
                k1 = max(1,k2-segLen);
                if k2~=1 && sum(isnan(orient0(k1:k2)))==0
                    welch0 = pwelch(orient0(k1:k2),[],[],freqs,50);
                else
                end
            elseif t2s(iM,iS)==-20 % missed spot
                seg = find(dnT0<100);
                if length(seg)>segLen
                    dseg = diff(seg);
                    breakit = cat(1,1,find(dseg>50),length(seg));
                    % use first 5s segment within 50cm of spot
                    keep = diff(breakit)>=700;
                    if any(keep)
                        k12 = seg(breakit(find(keep==1,1,'first')):breakit(find(keep==1,1,'first')+1)); %#ok<FNDSB>
                        k2 = breakit(find(keep==1,1,'first')+1);
                        welch1 = pwelch(orient0(k12(1:segLen)),[],[],freqs,50);
                    else
                    end
                end
                
                if ~isnan(k2) % get 5s segment far from spot / after finding spot
                    k4 = k2:length(orient0);
                    seg = find(dnT0>=100);
                    if length(seg)>segLen
                        dseg = diff(seg);
                        breakit = cat(1,1,find(dseg>50),length(seg));
                        keep = diff(breakit)>=segLen;
                        if any(keep)
                            k12 = seg(breakit(find(keep==1,1,'first')):breakit(find(keep==1,1,'first')+1)); %#ok<FNDSB>
                            welch2 = pwelch(orient0(k12(1:segLen)),[],[],freqs,50);
                        else
                        end
                    else
                    end
                else
                end
            else
            end
            
            welchSpot = cat(1,welchSpot,welch0);
            welchCont = cat(1,welchCont,welch1);
            welchFar = cat(1,welchFar,welch2);
        end
    end
    disp(iM);
end


% mseb(freqs,nanmean(welchCont),nanstderr(welchCont));
% lineProps.col = {'r'}; mseb(freqs,nanmean(welchSpot),nanstderr(welchSpot),lineProps);
% lineProps.col = {'k'}; mseb(freqs,nanmean(welchFar),nanstderr(welchFar),lineProps);
plot(freqs,nanmean(log10(welchCont)),'color','b'); hold on;
plot(freqs,nanmean(log10(welchFar)),'color','k'); 
plot(freqs,nanmean(log10(welchSpot)),'color','r');
%figure(1); clf;
%lineProps.col = {'b'};
%mseb(freqs,nanmean(welchSpot(spotfound==0,:),1)',nanstderr(welchSpot(spotfound==0,:),[],1),lineProps);
%lineProps.col = {'r'};
%mseb(freqs,nanmean(welchSpot(spotfound==1,:),1)',nanstderr(welchSpot(spotfound==1,:),[],1),lineProps);
%lineProps.col = {'k'};
%mseb(freqs,nanmean(welchCont(spotfound==0,:),1)',nanstderr(welchCont(spotfound==0,:),[],1),lineProps);
%lineProps.col = {'g'};
%mseb(freqs,nanmean(welchCont(spotfound==1,:),1)',nanstderr(welchCont(spotfound==1,:),[],1),lineProps);
% plot(freqs,nanmean(log10(welchSpot(spotfound==0,:)),1),'b','linewidth',2); hold on;
% plot(freqs,nanmean(log10(welchSpot(spotfound==1,:)),1),'r','linewidth',2);
% plot(freqs,nanmean(log10(welchCont(spotfound==0,:)),1),'k','linewidth',2);
% plot(freqs,nanmean(log10(welchCont(spotfound==1,:)),1),'g','linewidth',2);
% legend('spot missed','spot found','control - missed','control - found');
