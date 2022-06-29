%parse2

condition = {'beforespot','afterspot','notfound','conc=0.1','conc=1 | 2'};
cond = 1;

spectype = {'Dline','dphi'};
spect = 1;

nM = max(mouse);
nS = max(sess);


ndb = 50;
ntb = 50;
nvb = 50;
dbin = quantile(dnT(~edge),ndb);
tbin = quantile(frame(~edge)/50,ntb);
vbin = quantile(nV(~edge),nvb);
% dbin = linspace(0,120,ndb);
% tbin = linspace(0,1500,ntb);
% vbin = linspace(0,40,nvb);
nF = length(freq);
switch spect
    case 1
        zF = nan(size(S));
    case 2
        zF = nan(size(S2));
end
for iF=1:nF
    switch spect
        case 1
            zF(iF,:)=nanzscore(abs(S(iF,:)));
        case 2
            zF(iF,:)=nanzscore(abs(S2(iF,:)));
    end
end

Stime = nan(length(freqset),ntb,nM,nS);
Sdist = nan(length(freqset),ndb,nM,nS);
Svel = nan(length(freqset),nvb,nM,nS);
nVD = [];
dnTD = [];
tsD = [];

for iM=1:nM
    for iS=1:nS
        if ALorKP0(iM,iS)==0
            ksess = 1;
            %         switch cond
            %             case 1
            %                 ksess = t2sF(iM,iS)>1;
            %             case 2
            %                 ksess = t2s(iM,iS)>1;
            %             case 3
            %                 ksess = t2s(iM,iS)==-20;
            %             case 4
            %                 ksess = t2s(iM,iS)>1 & conc0(iM,iS)==0;
            %             case 5
            %                 ksess = t2s(iM,iS)>1 & conc0(iM,iS)>0;
            %         end
            if ksess
                k = mouse==iM & sess==iS;
                switch spect
                    case 1
                        kS = mouse1==iM & sess1==iS;
                    case 2
                        kS = mouse2==iM & sess2==iS;
                end
                switch cond
                    case 1
                        toss = edge(k);
                    case 2
                        toss = edge(k) | isstopped(k) | t2s0(k)==1;
                    case 3
                        toss = edge(k) | isstopped(k) | spotfound0trunc(k)==0;
                    case 4
                        toss = edge(k) | isstopped(k) | t2s0(k)==0;
                    case 5
                        toss = edge(k) | isstopped(k) | t2s0(k)==0;
                end
                dnTD0 = dnT(k);
                nVD0 = nV(k);
                tsD0 = frame(k)/50;
                
                tsD0(toss) = [];
                dnTD0(toss) = [];
                nVD0(toss) = [];
                
                switch spect
                    case 1
                        tsS0 = times(kS);
                    case 2
                        tsS0 = times2(kS);
                end
                if length(tsS0) > 1
                    
                    tsS1 = unique(interp1(tsS0,tsS0,tsD0,'nearest'));
                    
                    idx = unique(interp1(tsS0,1:length(tsS0),tsS1,'nearest'));
                    tsS1(idx==0 | isnan(idx))=[];
                    idx(idx==0 | isnan(idx))=[];
                    
                    if ~isempty(idx) && length(tsS1)>1
                        
                        S0 = zF(:,kS);
                        
                        %                     switch cond
                        %                         case {1,4,5}
                        %                             [~,~,tH] = histcounts(t2s(iM,iS)+1-tsS1,tbin); % +1s beyond capture
                        %                         case 2
                        %                             [~,~,tH] = histcounts(tsS1+1+t2s(iM,iS),tbin); % after capture
                        %                         case 3
                        %                             t2match = 2.5*nanmedian(t2s(iM,:)>1);
                        %                             tmatch = find(tsD0>t2match,1,'first');
                        %                             [~,~,tH] = histcounts(tmatch-tsS1+1,tbin); % get matched time to last found spot
                        %                     end
                        %
                        %                     tsST = tsS1;
                        %                     tsST(tH==0)=[];
                        %
                        %                     nT = length(tbin);
                        %                     tix = 1:nT;
                        %                     for iT=1:nT
                        %                         idxT = unique(interp1(idx,idx,find(tH==iT),'nearest'));
                        %                         tsST(isnan(idxT))=[];
                        %                         idxT(isnan(idxT))=[];
                        %                         if ~isempty(idxT)
                        %                             Stime(:,tix(iT),iM,iS)=nanmean(S0(:,idxT),2);
                        %                             tsD = cat(1,tsD,tsST);
                        %                         end
                        %                     end
                        
                        
                        dnTS = unique(interp1(tsD0,dnTD0,tsS1,'nearest'));
                        [~,~,dH] = histcounts(dnTS,dbin);
                        dnTS(dH==0)=[];
                        dH(dH==0)=[];
                        
                        nT = length(dbin);
                        dix = 1:nT;
                        for iT=1:nT
                            idxD = unique(interp1(idx,idx,find(dH==iT),'nearest'));
                            dnTS(isnan(idxD))=[];
                            idxD(isnan(idxD))=[];
                            if ~isempty(idxD)
                                Sdist(:,dix(iT),iM,iS)=nanmean(S0(:,idxD),2);
                                dnTD = cat(1,dnTD,dnTS);
                            end
                        end
                        
                        
                        nVS = unique(interp1(tsD0,nVD0,tsS1,'nearest'));
                        [~,~,bnV0] = histcounts(nVS,vbin);
                        nVS(bnV0==0)=[];
                        bnV0(bnV0==0)=[];
                        
                        nT = length(vbin);
                        vix = 1:nT;
                        for iT=1:nT
                            idxV = interp1(idx,idx,find(bnV0==iT),'nearest');
                            nVS(isnan(idxV))=[];
                            idxV(isnan(idxV))=[];
                            if ~isempty(idxV)
                                Svel(:,vix(iT),iM,iS)=nanmean(S0(:,idxV),2);
                                nVD = cat(1,nVD,nVS);
                            end
                        end
                    end
                end
            end
        end
    end
    disp(iM);
end


dS1 = nanmean(nanmean(Sdist,4),3);
tS1 = nanmean(nanmean(Stime,4),3);
vS1 = nanmean(nanmean(Svel,4),3);

figure;
pcolor(tbin,freq,smooth2a(tS1,2,1)); shading flat; xlabel('t-tspot'); ylabel('freq (Hz)'); title('PSD'); colorbar;
figure;
pcolor(dbin,freq,smooth2a(dS1,2,1)); shading flat; xlabel('dist (cm)'); ylabel('freq (Hz)'); title('PSD'); colorbar;
figure;
pcolor(vbin,freq,smooth2a(vS1,2,1)); shading flat; xlabel('v (cm/s)'); ylabel('freq (Hz)'); title('PSD'); colorbar;





