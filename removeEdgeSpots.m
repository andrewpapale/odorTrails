% removeEdgeSpots
% 2018-09-17 AndyP

edgespot = zeros(size(edge));
cxr = ones(size(cx0));
cyr = ones(size(cy0));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if t2s(iM,iS)==-40
            cxr(iM,iS)=0;
            cyr(iM,iS)=0;
            edgespot(k)=1;
        end
    end
end
            