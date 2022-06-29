% stim PETH
PETHx = nan(size(x));
PETHt = nan(size(x));
PETHs = nan(size(x));
dF = 200;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kL = mouseL ==iM & sessL ==iS;
        if sum(k)>0
            tStart = find(k==1,1,'first');
            sF0 = find(stimF(k)==1);
            sF0(sF0 < dF+1)=[];
            sF0(sF0 > find(k==1,1,'last')-dF+1)=[];
            %k0 = cat(1,1,diff(sF0)>dF+1);
            %sF0(~k0)=[];
            a0 = ampOut(kL);
            for iF=1:length(sF0)
                t0 = t(tStart+(sF0(iF)-dF-1):tStart+(sF0(iF)+dF-1))-t(sF0(iF)-1); % time from stim
                dvec = -dF:dF;
                if any(diff(t0)<0)
                    [~,trt]=max(t0);
                    dvec = dvec(trt:end);
                    t0 = t0(trt:end);
                else
                    trt = 1;
                end
                %[~,bt]=histcounts(t0,linspace(t0(1),t0(end),50));
                PETHx(trt+tStart+(sF0(iF)-dF-2):tStart+(sF0(iF)+dF-1))=dvec; % frame from stim
                PETHt(trt+tStart+(sF0(iF)-dF-2):tStart+(sF0(iF)+dF-1))=t0; % put times into time matrix
                PETHs(trt+tStart+(sF0(iF)-dF-2):tStart+(sF0(iF)+dF-1))=a0(iF);
            end
        end
    end
end

