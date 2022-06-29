
% 2018-02-21 AndyP 
% compute a logistic regression

k = t2s~=-40 & (t2s>1 | t2s==-20);
[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf));
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
X = [ALorKP0(k),conc0(k),trial0(k),bait0(k),bS(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),meandphi(k),fractoward(k),angcast0(k)];
%X = [ALorKP0(k),bait0(k),conc0(k)];
%X = bait0(k); % yes
%X = ALorKP0(k); % no
%X = conc0(k); % no
%X = mouse0(k);% yes
%X = bS(k); % no
%X = initD(k); % maybe
%X = initO(k); % no
%X = meanV(k); % no
%X = meannV(k); % no
%X = meandphi(k); % no
%X = conc1(k);
%X = [bait0(k),meanV(k),initD(k)];
%rn = randi(10,size(k));
%X = [rn(k)];
%X = fractoward(k);
%X = pE(k);
Y = spotfound(k);

nShuf = 500;
k1 = ones(size(Y));
ytest = nan(nShuf,1);
ytrain = nan(nShuf,1);
for iShuf = 1:nShuf
    rng('shuffle');
    train = zeros(size(k1));
    ri = randi(length(k1),[round(0.5*length(k1)),1]);
    train(ri)=1;
    test = ~train;
    %disp(any(train & test))
    [b,~,~] = glmfit(X(train==1),[Y(train==1),ones(sum(train==1),1)],'binomial','link','logit');
    ytrain0 = nansum(Y(train==1))./sum(train==1);
    
    ytest0 = glmval(b(1:end),X(test==1),'logit','constant','on','size',ones(sum(test==1),1));
    
    ytest(iShuf)=nanmedian(ytest0);
    ytrain(iShuf)= ytrain0;
    disp(iShuf);
end

bins = linspace(0,1,30);
H1 = hist(ytrain,bins);
H2 = hist(ytest,bins);
bar(linspace(0,1,30),cat(1,H1,H2)');
%[h,p,ks2stat] = kstest2(ytrain,ytest)
[h,p,ci,stats] = ttest2(ytrain,ytest)


