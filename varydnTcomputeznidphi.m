% vary dnT, compute znidphi as f(n) of concentration

bins = linspace(0.1,10,10);

H = [];
for iT=bins
    k = dnT<= iT;
    H0 = histcn(conc(k==1),2:3,'AccumData',znidphi(k==1),'fun',@nanmean);
    H = cat(2,H,H0);
end