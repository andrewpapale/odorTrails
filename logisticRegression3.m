
% 2018-02-21 AndyP 
% compute a logistic regression

b = [];
p = [];
for iR=1:15

k = t2s~=-40 & ThrV==iR & ~isinf(meandphiV);
%X = [ALorKP0(k),conc0(k),trial0(k),bait0(k),bS(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),meandphi(k)];
%X = [ALorKP0(k),bait0(k),conc0(k)];
% X = bait0(k); % yes
% X = ALorKP0(k); % no
% X = conc0(k); % no
%X = mouse0(k);% yes
% X = bS(k); % no
% X = initD(k); % maybe
% X = initO(k); % no
%X = meanV(k); % no
%X = meannV(k); % no
% X = meandphi(k); % no
%X = conc1(k);
% X = [bait0(k),mouse0(k),initD(k)];
%rn = randi(10,size(k));
%X = [rn(k)];
X = [ALorKP0(k),bait0(k),conc0(k),trial0(k),sess0(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),meandphi(k),meandphiV(k),fractoward(k),fracedge(k),ddnV(k),ddV(k),spotring(k)];
Y = spotfound(k);
[b0,~,stats] = glmfit(X,Y,'binomial','link','logit','options',options);
b(:,iR)=b0;
p(:,iR)=stats.p;
end

% nShuf = 200;
% k1 = ones(size(Y));
% ytest = nan(nShuf,1);
% ytrain = nan(nShuf,1);
% for iShuf = 1:nShuf
%     train = zeros(size(k1));
%     ri = randi(length(k1),[round(0.5*length(k1)),1]);
%     train(ri)=1;
%     test = ~train;
%     %disp(any(train & test))
%     [b,~,~] = glmfit(X(train==1),[Y(train==1),ones(sum(train==1),1)],'binomial','link','logit');
%     ytrain0 = nansum(Y(train==1))./sum(train==1);
%     
%     ytest0 = glmval(b(2:end),X(test==1),'logit','constant','off','size',ones(sum(test==1),1));
%     
%     ytest(iShuf)=nanmean(ytest0);
%     ytrain(iShuf)= ytrain0;
%     disp(iShuf);
% end
% 
% bins = linspace(0,1,30);
% H1 = hist(ytrain,bins);
% H2 = hist(ytest,bins);
% bar(linspace(0,1,30),cat(1,H1,H2)');
% %[h,p,ks2stat] = kstest2(ytrain,ytest)
% [h,p,ci,stats] = ttest2(ytrain,ytest)
F1 = figure(1); clf; %subplot(1,2,1);
plot(thrVec,log10(p(2:7,:)'),'linewidth',2);
xlabel('threshold distance (cm)');
ylabel('log10(p-val)');
legend('ALorKP','bait','conc','trial','session','mouse');
set(gca,'YLim',[-9,0]);

F2 = figure(2); clf; %subplot(1,2,1);
plot(thrVec,log10(p(8:13,:)'),'linewidth',2);
xlabel('threshold distance (cm)');
ylabel('log10(p-val)');
legend('initD','initO','meanV','meannV','meandphi','meandphi/V');
set(gca,'YLim',[-9,0]);

F3 = figure(3); clf; % subplot(1,2,1);
plot(thrVec,log10(p(14:18,:)'),'linewidth',2);
xlabel('threshold distance (cm)');
ylabel('log10(p-val)');
legend('fractoward','fracedge','ddnV','ddV','spotring');
set(gca,'YLim',[-9,0]);