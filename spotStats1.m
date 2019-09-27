% 2018-09-14 AndyP
% Figure2A 
% logistic regression on % correct
k = t2s>mint2s | t2s==-20;
[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf)); % bin session
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);

X = [ALorKP0(k),conc0(k),trial0(k),bait0(k),bS(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),meandphi(k),fractoward(k),angcast0(k)];



