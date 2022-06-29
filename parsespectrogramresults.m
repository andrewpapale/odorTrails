% parsespectrogram results

nM = max(mouse);
nS = max(sess);


ndb = 30;
ntb = 30;
nvb = 30;
dbin = quantile(dT(t2s0==1 & ~edge & ~isstopped),ndb);
tbin = quantile(frame(t2s0==1 & ~edge & ~isstopped)/50,ntb);
vbin = quantile(nV(t2s0==1 & ~edge & ~isstopped),nvb);
%%
zS = nan(size(S));
zS2 = nan(size(S2));
%for iM=1:nM
    %for iS=1:nS
    %k = mouse1==iM;
    k2 = mouse2==iM;
       % k = mouse1==iM & sess1==iS;
       % k2 = mouse3==iM & sess3==iS;
        for iF=1:size(S,1)
            kF = S(iF,:);
            %kF2 = S2(iF,k2);
            zS0 = nanzscore(abs(kF));
            %zS2 = nanzscore(abs(kF2));
            zS(iF,:)=zS0;
            %zS2(iF,k2)=zS2;
        end
    %end
    disp(iM);
%end


%%
Stime = nan(100,ntb,nM,nS);
Sdist = nan(100,ndb,nM,nS);
vS1 = nan(100,nvb,nM,nS);
S1 = [];
nV1 = [];
mouse3 = [];
sess3 = [];
dnT1 = [];
frame2 = [];

for iM=1:nM
    for iS=1:nS
        
        k = mouse==iM & sess==iS;
        firstF = find(k==1,1,'first'); % index to first frame in the session iS,iM
        if sum(k)>0 && t2s(iM,iS)>1 % if successful session with data
            
            kS = mouse1==iM & sess1==iS; % find windows for iS, iM
            
            frame0 = 1:sum(k); % create vector for each frame, starting from 1
            toss = edge(k) | isstopped(k) | t2s0(k)==0;
            dnTx = dnT(k);
            nVx = nV(k);
            dnTx(toss) = nan;
            nVx(toss) = nan;
            
            frame1 = round(times(kS)*50); % convert spectrogram times -> frames
            nF = length(frame0);
            for iF=1:nF % for each data frame
                tempframe = frame1(frame1==frame0(iF)); % find spec frame==data frame iF 
                if ~isempty(tempframe)
                    frame2 = cat(1,frame2,tempframe); % save frame as frame2
                end
            end
            frame2 = frame2(~isnan(frame2)); 
            frame2 = unique(frame2); % get unique frames for binning
            frame2(frame2 > max(frame0))=[]; % ??? why do i need to do this
            nF = length(frame2); 
            indx = [];
            for iF=1:nF % for each unique matched frame
                indx0 = find(frame1==frame2(iF),1,'first'); % find index of spectrogram equivalent to matched frame
                if ~isempty(indx0)
                    indx = cat(1,indx,indx0);
                else
                    frame2(iF)=-1;
                end
            end
            frame2(frame2==-1)=[];
            
            S0 = S(:,kS);
            S0 = S0(:,indx);
            S1 = cat(2,S1,S0);
            % order by time bin relative to spot
            spotframe = round(t2s(iM,iS)*50)+0; % spotframe
            [~,tH] = histc((spotframe-frame1(indx))/50,tbin);
            [tH,six] = sort(tH,'ascend');
            ntH = unique(tH);
            nT = length(ntH);
            for iT=2:nT
                    Stime(:,ntH(iT),iM,iS)=nanmean(S0(:,six(tH==ntH(iT))),2);
            end
            
            % order by distance bin from spot
            % get closest non-nan point
            nV0 = nVx(frame2);
            nV1 = cat(1,nV1,nV0);
            dnT0 = dnTx(frame2);
            dnT1 = cat(1,dnT1,dnT0);
            mouse3 = cat(1,mouse3,repmat(iM,[length(frame2),1]));
            sess3 = cat(1,sess3,repmat(iS,[length(frame2),1]));
            [~,dH] = histc(dnT0,dbin);
            [dH,six] = sort(dH,'ascend');
            ndH = unique(dH);
            nT = length(ndH);
            for iT=2:nT
                %disp(sum(six(dH==ndH(iT))));
                    Sdist(:,ndH(iT),iM,iS)=nanmean(S0(:,six(dH==ndH(iT))),2);
            end
            
            [~,bnV0] = histc(nV0,vbin);
            [bnV0,six]=sort(bnV0,'ascend');
            nvH = unique(bnV0);
            for iT=2:length(unique(bnV0))
                vS1(1:100,nvH(iT),iM,iS)=nanmean(S1(:,six(bnV0==nvH(iT))),2);
            end
            
            
        end
    end
    disp(iM);
end
dS = nanmean(nanmean(Sdist,4),3);
tS = nanmean(nanmean(Stime,4),3);
vS2 = nanmean(nanmean(vS1,4),3);

