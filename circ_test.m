% M = 10:1:1000;
% pval = nan(length(M),1);
% F = nan(length(M),1);
% n = nan(length(M),1);
% rw = nan(length(M),1);
% for iB=M
%     bins = linspace(-pi,pi,iB);
%     [MO,bO]=histc(initO(:)*pi/180,bins);
%     b1 = nan(size(bO));
%     w1 = nan(size(bO));
%     k = bO>0;
%     w1(k)=MO(bO(k));
%     b1(k)=bins(bO(k));
%     [pval0, table0,n0,rw0] = circ_wwtest(cat(1,b1(k1),b1(k2)),cat(1,ones(sum(k1(:)),1),ones(sum(k2(:)),1)+1),cat(1,w1(k1),w1(k2)));
%     pval(iB,1)=pval0;
%     F(iB,1)=table{2,5};
%     n(iB,1)=n0;
%     rw(iB,1)=rw0;
% end


k1 = ALorKP0==0 & ~isnan(initO);
k2 = ALorKP0==1 & ~isnan(initO);
[pval0, k0, ~] = circ_kuipertest(initO(k1)*pi/180, initO(k2)*pi/180, 100, false);

nM = 300;
pval = nan(nM,1);
kval = nan(nM,1);
initO1 = initO(~isnan(initO));
for iB=1:nM
    rng('shuffle');
    idx = 1:length(initO1);
    k1 = randi(length(initO1),[300,1]);
    idx2 = [];
    for iL=1:length(k1)
        idx1 = find(idx==k1(iL),1,'first');
        idx2 = cat(1,idx2,idx1);
    end
    idx(idx2)=[];
    k2 = randi(length(idx),[300,1]);
    k2 = idx(k2);
    
    [pval0, k0, ~] = circ_kuipertest(initO1(k1)*pi/180, initO1(k2)*pi/180, 100, false);
    pval(iB,1)=pval0;
    kval(iB,1)=k0;
end