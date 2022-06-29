function x1 = getNearestNonNan(x,t,tol)

t2 = find(~isnan(x));

x1 = nan(size(t));
for iT=1:length(t)
    [minT,Ix] = min(abs(t2-t(iT)));
    if minT <= tol
        x1(iT) = x(t2(Ix));
    end
end

end