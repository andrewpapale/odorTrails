mouse2 = [];
welchSpot = [];
welchCont = [];
kSpotS = [];
kSpotF = [];
kStartC = [];
kStopC = [];
dntspot = [];
dntcont = [];
freqs = linspace(0.1,10,100);
spotfound = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS & ~edge;
        
        if sum(k)>0
            t2s0 = t2s(iM,iS);
            dnT0 = dnT(k);
            
            if t2s0>0
                t2s0 = round(t2s0*50);
                flag = 0;
                seg = [];
            elseif t2s0==-20
                seg = find(~isnan(dnT0(1:500)));
            end
            if ~isempty(seg) && t2s(iM,iS)>-40
                
                alignit = seg(1)-1;
                dseg = diff(seg);
                breakit = cat(1,1,find(dseg>1),length(seg));
                [maxseg,Imaxseg] = max(diff(breakit));
                
                distcont = find(dnT0<20,1,'first');
                if ~isempty(distcont)
                    t2s0 = distcont;
                    flag = 0;
                else
                    t2s0 = nan;
                    flag = 1;
                end
                
            elseif t2s0==-40
                flag = 1;
            end
            
            
            
            orient0 = orient(k);
            
            
            if t2s0 > 50 & ~flag %#ok<AND2>
                
                if t2s(iM,iS)>0
                    spotfound = cat(1,spotfound,1);
                elseif t2s(iM,iS)==-20
                    spotfound = cat(1,spotfound,0);
                else
                    keyboard;
                end
                
                
                k1 = min(t2s0,length(orient0));
                k2 = max(1,t2s0-500);
                kSpotS = cat(1,kSpotS,k1);
                kSpotF = cat(1,kSpotF,k2);
                dntspot = cat(1,dntspot,nanmean(dnT0(k2:k1)));
                welch0 =pwelch(orient0(k2:k1),[],[],freqs,50);
                welchSpot = cat(1,welchSpot,welch0);
                mouse2 = cat(1,mouse2,repmat(iM,[1,length(freqs)]));
                
                
                % find equivalent length segment not going to spot
                
                seg = find(~isnan(orient0(t2s0+500:end)));
                
                if ~isempty(seg)
                    
                    alignit = seg(1)-1;
                    dseg = diff(seg);
                    breakit = cat(1,1,find(dseg>1),length(seg));
                    [maxseg,Imaxseg] = max(diff(breakit));
                    
                    if maxseg > 50
                        fstart = breakit(Imaxseg)+t2s0+500+alignit;
                        kStartC = cat(1,kStartC,fstart);
                        fstop = min(cat(1,length(orient0),breakit(Imaxseg+1)+t2s0+500+alignit,fstart+length(k2:k1)));
                        kStopC = cat(1,kStopC,fstop);
                        dntcont = cat(1,dntcont,nanmean(dnT0(fstart:fstop)));
                        orient1 = orient0(fstart:fstop);
                        welch1 = pwelch(orient1,[],[],freqs,50);
                    else
                        welch1 = nan(1,length(freqs));
                        dntcont = cat(1,dntcont,nan);
                        kStopC = cat(1,kStopC,nan);
                        kStartC = cat(1,kStartC,nan);
                    end
                    % keyboard;
                else
                    dntcont = cat(1,dntcont,nan);
                    kStopC = cat(1,kStopC,nan);
                    kStartC = cat(1,kStartC,nan);
                    welch1 = nan(1,length(freqs));
                end
                welchCont = cat(1,welchCont,welch1);
            else
                
            end
        end
    end
    disp(iM);
end


figure(1); clf;
%lineProps.col = {'b'};
%mseb(freqs,nanmean(welchSpot(spotfound==0,:),1)',nanstderr(welchSpot(spotfound==0,:),[],1),lineProps);
%lineProps.col = {'r'};
%mseb(freqs,nanmean(welchSpot(spotfound==1,:),1)',nanstderr(welchSpot(spotfound==1,:),[],1),lineProps);
%lineProps.col = {'k'};
%mseb(freqs,nanmean(welchCont(spotfound==0,:),1)',nanstderr(welchCont(spotfound==0,:),[],1),lineProps);
%lineProps.col = {'g'};
%mseb(freqs,nanmean(welchCont(spotfound==1,:),1)',nanstderr(welchCont(spotfound==1,:),[],1),lineProps);
plot(freqs,nanmean(log10(welchSpot(spotfound==0,:)),1),'b','linewidth',2); hold on;
plot(freqs,nanmean(log10(welchSpot(spotfound==1,:)),1),'r','linewidth',2);
plot(freqs,nanmean(log10(welchCont(spotfound==0,:)),1),'k','linewidth',2);
plot(freqs,nanmean(log10(welchCont(spotfound==1,:)),1),'g','linewidth',2);
legend('spot missed','spot found','control - missed','control - found');
