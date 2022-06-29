function tdist = dist2nan(y)
% compute distance to nearest nan

tnan = find(isnan(y));
tdist = nan(size(y));
for iT=1:length(tdist)
    tdist(iT)=min(abs(iT-tnan));
end
