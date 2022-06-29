
% 2018-02-21 AndyP 
% compute a logistic regression

b = [];
p = [];
F = figure(1); clf;
for iR=1:15

k = t2s~=-40 & ThrV==iR & ~isinf(meandphiV);
X = [ALorKP0(k),conc0(k),trial0(k),bait0(k),bS(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),meandphi(k),pE(k)];
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
X = [ALorKP0(k),bait0(k),fractoward(k)];
Y = spotfound(k);
% [b0,~,stats] = glmfit(X,Y,'binomial','link','logit','options',options);
% b(:,iR)=b0;
% p(:,iR)=stats.p;
mdl = fitglm(X,Y,'Distribution','binomial','Link','logit');
score_log = mdl.Fitted.Probability;
[Xlog,Ylog,Tlog,AUClog] = perfcurve(Y,score_log,'1');
switch iR
    case {1,2,3,4,5}
        subplot(1,3,1); plot(Xlog,Ylog); hold on;
    case {6,7,8,9,10}
        subplot(1,3,2); plot(Xlog,Ylog); hold on;
    case {11,12,13,14,15}
        subplot(1,3,3); plot(Xlog,Ylog); hold on;
end
end
subplot(1,3,1); legend(mat2str(thrVec(1)),mat2str(thrVec(2)),mat2str(thrVec(3)),mat2str(thrVec(4)),mat2str(thrVec(5)));
subplot(1,3,2); legend(mat2str(thrVec(6)),mat2str(thrVec(7)),mat2str(thrVec(8)),mat2str(thrVec(9)),mat2str(thrVec(10)));
subplot(1,3,3); legend(mat2str(thrVec(11)),mat2str(thrVec(12)),mat2str(thrVec(13)),mat2str(thrVec(14)),mat2str(thrVec(15)));

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
% F1 = figure(1); subplot(1,2,1);
% plot(thrVec,log10(p(2:5,:)'),'linewidth',2);
% xlabel('threshold distance (cm)');
% ylabel('log10(p-val)');
% legend('conc','trial','session','mouse');
% set(gca,'YLim',[-7,0]);
% 
% F2 = figure(2); subplot(1,2,1);
% plot(thrVec,log10(p(6:11,:)'),'linewidth',2);
% xlabel('threshold distance (cm)');
% ylabel('log10(p-val)');
% legend('initD','initO','meanV','meannV','meandphi','meandphi/V');
% set(gca,'YLim',[-7,0]);
% 
% F3 = figure(3); subplot(1,2,1);
% plot(thrVec,log10(p(12:16,:)'),'linewidth',2);
% xlabel('threshold distance (cm)');
% ylabel('log10(p-val)');
% legend('fractoward','fracedge','ddnV','ddV','spotring');
% set(gca,'YLim',[-7,0]);