ftrail = dir('*-Odor-Trail.mat');
fcp = dir('*-ConnectingPoint.mat');

assert(length(ftrail)==length(fcp),'file mismatch');

nF = length(fcp);

for iF=1:nF
    
    load(ftrail(iF).name);
    load(fcp(iF).name);
    
    saveStr = fcp(iF).name;
    disp(saveStr);
    save(saveStr,'xC','yC','outside','inside','vec','xC0','yC0','arm');
end