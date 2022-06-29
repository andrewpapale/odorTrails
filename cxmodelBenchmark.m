%cxmodelbenchmark
%
%Requires Matlab parallel computing toolbox, cxmodel5ps0.m and
%network05-0.mat
%
%Execution should take about 30-180 minutes.
%
%R.V. Williams-Garcia 2017

[~,hostName]= system('hostname');

p = gcp;
if isempty(p)
    poolSize = num2str(0);
else
    poolSize = num2str(p.NumWorkers);
end
delete(gcp)

rngState = uint32(1583493145);
rng(rngState)

kappa = 1.3;
load('network05-0.mat','sig_del','weightmat')
A = kappa*weightmat;
A(A>1) = 1;
N = length(weightmat);
tic
[cShapes,~,~] = cxmodel5ps0(A,1E5,ones(N,1),sig_del,zeros(N),1E5);
runtime = toc;

cSizes = cellfun(@(x) size(x,1),cShapes);
cDurations = cellfun(@(x) x(end,2)-x(1,2)+1,cShapes);
        
save(strcat('Benchmark',hostName,poolSize),'cDurations','cSizes','kappa','rngState','runtime')
