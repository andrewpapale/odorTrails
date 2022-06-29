cnx = nan(size(x));
cny = nan(size(y));
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k) > 0
            cnx(k) = nx(k)-cx1(k);
            cny(k) = nx(k)-cy1(k);
        end
    end
end
