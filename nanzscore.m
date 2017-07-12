function y = nanzscore(x)

y = (x - nanmean(x))./nanstd(x);


end

