function zmC = zScoreMouse(zC,mouse)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
nM = max(mouse);

zmC = nan(size(zC));
for iM=1:nM
    k = mouse==iM;
    if nansum(k)>0
        zmC(k) = nanzscore(zC(k));
    end
end


end

