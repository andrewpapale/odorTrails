function log10nC = logVar(nC)

min = 1E-10;
max = 1E10;

nC(nC<=min)=nan;
nC(nC>=max)=nan;
log10nC = log10(nC);

end