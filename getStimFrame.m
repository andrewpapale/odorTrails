% getStimFrame
nM = 5;
nS = 50;
stimF = zeros(size(x));
f0 = nan(size(x));
%ISI0 = nan(size(x));
%AMP0 = nan(size(x));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kL = mouseL ==iM & sessL ==iS;
        if sum(k)>0 && sum(kL)>0
            %ISI0(k)=repmat(ISI(iM,iS),[sum(k),1]);
            %AMP0(k)=repmat(Vout(iM,iS),[sum(k),1]);
            frame0 = 1:sum(k);
            firstT = find(k==1,1,'first');
            %sF0 = interp1(t(k),frame0,Stim(kL),'nearest'); %2019-07-19 AndyP from Annie's dataset
            sF0 = interp1(t(k),frame0,timOut(kL),'nearest');
            stimF(firstT+sF0)=1; % stim occurs ~1 frame after timestamp
            f0(firstT+sF0)=ampOut(kL);
        end
    end
end