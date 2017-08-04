function log10nC = logVar(nC)

min = 1E-10;
max = 1E10;

nC(nC<=min)=min;
nC(nC>=max)=max;
log10nC = log10(nC);

end