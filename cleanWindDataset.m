% clean wind dataset
% 2019-02-28 AndyP

% remove spots or fans in the wrong position
% after running in WindZone...

remove = zeros(size(x));
removeT = zeros(size(fff));
removeS = zeros(nM,nS);


% get wP
wP0 = nan(size(x));
wP = nan(nM,nS);
for iM=1:nM
    for iS=1:nS
        fx1 = fx0(iM,iS);
        fy1 = fy0(iM,iS);
        k = mouse==iM & sess==iS;
        if sum(k)>0
            if fy1 > 250 & fx1 < 450 & fx1 > 300
                wP(iM,iS)=1;
            elseif fy1 > 250 & fx1 > 450
                wP(iM,iS)=2;
            elseif fy1 < 250 & fx1 > 450
                wP(iM,iS)=3;
            elseif fy1 < 250 & fx1 < 450 & fx1 > 300
                wP(iM,iS)=4;
            elseif fy1 < 250 & fx1 < 300 & fx1 > 125
                wP(iM,iS)=5;
            elseif fy1 < 250 & fx1 < 125
                wP(iM,iS)=6;
            elseif fy1 > 250 & fx1 < 125
                wP(iM,iS)=7;
            elseif fy1 > 250 & fx1 > 125 & fx1 < 300 %#ok<*AND2>
                wP(iM,iS)=8;
            else
            end
            wP0(k) = repmat(wP(iM,iS),[sum(k),1]);
        end
    end
end






ktoss = (wP0==1 & wind0==1) | (wP0==4 & wind0==1) | (wP0==7 & wind0==0) | (cx1 < 100 & cy1 < 200) | (cx1 < 140 & cy1 > 360) | (cx1 > 300 & cy1 > 380) | ...
    (mouse==2 & sess==12) | (mouse==2 & sess==13) | (mouse==2 & sess==14) | (mouse==2 & sess==44) | (mouse==2 & sess==73)| ...
    (mouse==4 & sess==5) | (cx1 > 300 & wP0==1);
remove(ktoss)=1;

s0 = unique(sess(ktoss));
m0 = unique(mouse(ktoss));

for iM=1:length(m0)
    for iS=1:length(s0)
        removeS(m0(iM),s0(iS)) = 1;
        ix = find(mouseT==iM & sessT==iS);
        removeT(ix)=1;
    end
end



