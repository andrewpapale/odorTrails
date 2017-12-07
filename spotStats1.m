% spot spot stats
% run after getTime2Spot and timeInRing

%
%
%

% time to spot (first time)

k = t2s > 0 & ~isnan(mouse0);

[p,table,stats]=anovan(t2s(k),{mouse0(k),conc0(k),bait0(k),ALorKP0(k)},'varnames',{'mouse','concentration','bait','ALorKP'},'model','interaction');

% strong effect of mouse => z-score within mouse
% interaction of mouse and baited/unbaited

zt2s = nan(size(t2s));
nR = size(mouse0,1);
for iR=1:nR
    k = t2s(iR,:)>0;
    zt2s(iR,k)=nanzscore(t2s(iR,k));
end

k = t2s>0 & ~isnan(mouse0);
[p,table,stats]=anovan(zt2s(k),{mouse0(k),conc0(k),bait0(k),ALorKP0(k)},'varnames',{'mouse','concentration','bait','ALorKP'},'model','interaction');
% no effects

% remove mouse, add session
[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf));
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
k = t2s > 0 & bS > 0 & ~isnan(mouse0);
[p,table,stats]=anovan(zt2s(k),{bS(k),conc0(k),bait0(k)},'varnames',{'session (binned)','concentration','bait/unbaited'},'model','interaction');
% no effect of session
% no effect of bait/unbait
boxplot(zt2s(k),bait0(k),'notch','on');
[p,h,stats] = ranksum(zt2s(k & bait0==0),zt2s(k & bait0==1));

boxplot(zt2s(k),{mouse0(k),bait0(k)},'notch','on');

k = zt2s > 0 & ~isnan(zt2s);
[p,table,stats]=anovan(zt2s(k),{conc0(k),bait0(k),ALorKP0(k),mouse0(k)},'varnames',{'concentration','bait/unbaited','ALorKP','mouse'},'model','interaction');

% initD
[p,table,stats]=anovan(t2s(k),{mouse0(k),conc0(k),bait0(k),ALorKP0(k),initD(k)},'varnames',{'mouse','concentration','bait','ALorKP','initD'},'model','interaction','continuous',5);

[p,table,stats]=anovan(zt2s(k)./initD(k),{mouse0(k),conc0(k),bait0(k),ALorKP0(k)},'varnames',{'mouse','concentration','bait','ALorKP'},'model','interaction');

%
%
%

% success / fail
hitmiss = nan(size(t2s));
hitmiss(t2s>0)=1;
hitmiss(t2s==-20)=2;  % for mnrfit needs to be > 0
% k = ~isnan(hitmiss);
% [tbl,chi2,p,labels] = crosstab(hitmiss(k),bait0(k));
% [h,p,stats] = fishertest(tbl);
% % mice are more successful when spot is baited
% 
% [tbl,chi2,p,labels] = crosstab(hitmiss(k),conc0(k));
% 
% % no effect of concentration
% 
% [tbl,chi2,p,labels] = crosstab(hitmiss(k),mouse0(k));
% % no effect of mouse
% 
% k = ~isnan(hitmiss) & bS>0;
% [tbl,chi2,p,labels] = crosstab(hitmiss(k),bS(k));
% effect of session ???

PrevTrial;

k = zt2s > 0 & ~isnan(zt2s);
[p,table,stats]=anovan(zt2s(k),{conc0(k),bait0(k),mouse0(k),prevbait(k),prevconc(k)},'varnames',{'concentration','bait/unbaited','mouse','prevbait','prevconc'},'model','interaction');


% multinomial logistic regression
k = ~isnan(hitmiss) & bS>0 & ~isnan(mouse0) & t2s > -40;
X = [bS(k),trial0(k),ALorKP0(k),bait0(k),conc0(k),mouse0(k),quadrant(k),prevbait(k),prevconc(k),initD(k),initO(k)];
Y = (hitmiss(k));
[B,dev,stats] = mnrfit(X,Y,'model','nominal');
[B,dev,stats] = mnrfit(X,Y,'model','ordinal'); % nominal or ordinal?
% ALorKP, bait, conc (I believe this is a false positive), and mouse are significant 


% success / fail
hitmiss = nan(size(t2s));
hitmiss(t2s>0)=1;
hitmiss(t2s==-20)=0;
% ALorKP
nansum(hitmiss(k) & ALorKP0(k)==0)./nansum(ALorKP0(k)==0)  % 80%
nansum(hitmiss(k) & ALorKP0(k)==1)./nansum(ALorKP0(k)==1)  % 50%
% bait
nansum(hitmiss(k) & bait0(k)==0)./nansum(bait0(k)==0)  % 45%
nansum(hitmiss(k) & bait0(k)==1)./nansum(bait0(k)==1)  % 77%
% mouse
nansum(hitmiss(k) & mouse0(k)==1)./nansum(mouse0(k)==1)  % 62
nansum(hitmiss(k) & mouse0(k)==2)./nansum(mouse0(k)==2)  % 60
nansum(hitmiss(k) & mouse0(k)==3)./nansum(mouse0(k)==3)  % 51
nansum(hitmiss(k) & mouse0(k)==4)./nansum(mouse0(k)==4)  % 67
nansum(hitmiss(k) & mouse0(k)==5)./nansum(mouse0(k)==5)  % 77

%mouse x ALorKP
% AL
nansum(hitmiss(k) & mouse0(k)==1 & ALorKP0(k)==0)./nansum(mouse0(k)==1 & ALorKP0(k)==0)  % 81
nansum(hitmiss(k) & mouse0(k)==2 & ALorKP0(k)==0)./nansum(mouse0(k)==2 & ALorKP0(k)==0)  % 78
nansum(hitmiss(k) & mouse0(k)==3 & ALorKP0(k)==0)./nansum(mouse0(k)==3 & ALorKP0(k)==0)  % 69 
nansum(hitmiss(k) & mouse0(k)==4 & ALorKP0(k)==0)./nansum(mouse0(k)==4 & ALorKP0(k)==0)  % 81
nansum(hitmiss(k) & mouse0(k)==5 & ALorKP0(k)==0)./nansum(mouse0(k)==5 & ALorKP0(k)==0)  % 87
% KP
nansum(hitmiss(k) & mouse0(k)==1 & ALorKP0(k)==1)./nansum(mouse0(k)==1 & ALorKP0(k)==1)  % 44
nansum(hitmiss(k) & mouse0(k)==2 & ALorKP0(k)==1)./nansum(mouse0(k)==2 & ALorKP0(k)==1)  % 51
nansum(hitmiss(k) & mouse0(k)==3 & ALorKP0(k)==1)./nansum(mouse0(k)==3 & ALorKP0(k)==1)  % 39 
nansum(hitmiss(k) & mouse0(k)==4 & ALorKP0(k)==1)./nansum(mouse0(k)==4 & ALorKP0(k)==1)  % 56
nansum(hitmiss(k) & mouse0(k)==5 & ALorKP0(k)==1)./nansum(mouse0(k)==5 & ALorKP0(k)==1)  % 63
% conclusion all mice have reduced performance for KP
% question, is this because there are more unbaited trials for KP?

nansum(bait0(k)==0 & ALorKP0(k)==0) % 152
nansum(bait0(k)==0 & ALorKP0(k)==1) % 199

% maybe, look at baited/unbaited separately
nansum(hitmiss(k) & bait0(k)==0 & ALorKP0(k)==0)./nansum(bait0(k)==0 & ALorKP0(k)==0) % 59
nansum(hitmiss(k) & bait0(k)==1 & ALorKP0(k)==0)./nansum(bait0(k)==1 & ALorKP0(k)==0) % 93

nansum(hitmiss(k) & bait0(k)==0 & ALorKP0(k)==1)./nansum(bait0(k)==0 & ALorKP0(k)==1) % 34
nansum(hitmiss(k) & bait0(k)==1 & ALorKP0(k)==1)./nansum(bait0(k)==1 & ALorKP0(k)==1) % 62

% conclusion: performance is reduced for both baited and unbaited
% conditions



%
%
%
% re-visiting spot

k = ~isnan(mouse0);
revisit = t2s2 > 0;
clf;
H = histcn([sess0(k),bait0(k)],cat(2,1,quantile(sess0(k),5),Inf),0:1,'AccumData',+(revisit),'fun',@nanmean);
dH = histcn([sess0(k),bait0(k)],cat(2,1,quantile(sess0(k),5),Inf),0:1,'AccumData',+(revisit),'fun',@nanstderr);
lineProps.col = {'b'};
mseb(cat(2,1,quantile(sess0(k),5)),H(:,1)',dH(:,1)',lineProps);
lineProps.col = {'r'};
mseb(cat(2,1,quantile(sess0(k),5)),H(:,2)',dH(:,2)',lineProps);
xlabel('session (binned)','fontsize',18);
ylabel('average number of spot returns','fontsize',18)
legend('unbaited','baited');

[tbl,chi2,p,labels] = crosstab(revisit(k),bait0(k));
[h,p,stats] = fishertest(tbl);
% mice revisit more when spots are baited

k = ~isnan(mouse0) & bait0==1;
[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf));
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
[tbl,chi2,p,labels] = crosstab(revisit(k),bS(k));

k = ~isnan(mouse0) & bait0==0;
[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf));
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
[tbl,chi2,p,labels] = crosstab(revisit(k),bS(k));
% mice revisit more previously baited spots on early sessions


revisit = +revisit;
revisit(revisit==0)=2;
% multinomial logistic regression
k = ~isnan(revisit) & bS>0 & ~isnan(mouse0) & t2s > -40;
X = [bS(k),trial0(k),ALorKP0(k),bait0(k),conc0(k),mouse0(k),quadrant(k),prevbait(k),prevconc(k)];
Y = (revisit(k));
[B,dev,stats] = mnrfit(X,Y,'model','ordinal');
% session (binned), trial, bait, and concentration are significant

% session
revisit(revisit==2)=0;
nansum(revisit(k) & bS(k)==1)./nansum(bS(k)==1) % 45
nansum(revisit(k) & bS(k)==2)./nansum(bS(k)==2) % 24
nansum(revisit(k) & bS(k)==3)./nansum(bS(k)==3) % 12
nansum(revisit(k) & bS(k)==4)./nansum(bS(k)==4) % 12
nansum(revisit(k) & bS(k)==5)./nansum(bS(k)==5) % 06

% % trial (I believe false positive, no clear trend, could be due to imbalances in bait/unbaited across trials)
% nansum(revisit(k) & trial0(k)==1)./nansum(trial0(k)==1) % 19
% nansum(revisit(k) & trial0(k)==2)./nansum(trial0(k)==2) % 13
% nansum(revisit(k) & trial0(k)==3)./nansum(trial0(k)==3) % 28
% nansum(revisit(k) & trial0(k)==4)./nansum(trial0(k)==4) % 10
% nansum(revisit(k) & trial0(k)==5)./nansum(trial0(k)==5) % 08

% bait/unbaited
nansum(revisit(k) & bait0(k)==0)./nansum(bait0(k)==0) % 08
nansum(revisit(k) & bait0(k)==1)./nansum(bait0(k)==1) % 24

% concentration (possibly more revisiting during 2%)
nansum(revisit(k) & conc0(k)==0)./nansum(conc0(k)==0) % 13
nansum(revisit(k) & conc0(k)==1)./nansum(conc0(k)==1) % 14
nansum(revisit(k) & conc0(k)==2)./nansum(conc0(k)==2) % 24

%
%
%
% dphi as a f'n of distance from spot



log10dphi = logVar(abs(nidphi)./nV);
znidphi1 = zScoreMouse(log10dphi,mouse);

edge = ~(nearX>50 & nearY>50 & nearX<nanmax(nearX)-50 & nearY<nanmax(nearY)-50);
k = nV>0.1 & ~edge;
bdT = cat(2,0,quantile(dnT(k),30),Inf);

[~,dT0]=histc(dnT,bdT);
[p,table,stats]=anova1(znidphi1(k),dT0(k));

clear p h stats
bootstat = nan(length(bdT),1);
dbootstat = nan(length(bdT),1);
for iD=1:length(bdT)-1
    disp(iD);
    k0 = k & (dnT >= bdT(iD) & dnT < bdT(iD+1));
    randznidphi = znidphi1(randi(length(znidphi1),[sum(k0),1]));
    bootstat(iD) = nanmedian(randznidphi);
    dbootstat(iD)= nanstderr(randznidphi);
    [p(iD),h(iD)] = ranksum(znidphi1(k0),randznidphi);
end

H = histcn([dnT(k)],bdT,'AccumData',znidphi1(k),'fun',@nanmedian);
dH = histcn([dnT(k)],bdT,'AccumData',znidphi1(k),'fun',@nanstderr);
N = histcn([dnT(k)],bdT);
clf;
lineProps.col = {'b'};
mseb(bdT(1:end-1),H(:,1)',dH(:,1)',lineProps);
line([0 90],[0 0],'color','k');

lineProps.col = {[0.5 0.5 0.5]};
mseb(bdT(1:end),bootstat',dbootstat',lineProps);

expfalsePos = length(bdT)-1;

for iD=1:length(bdT)-1
    if p(iD)<0.05/expfalsePos;
       text(bdT(iD),H(iD)+0.1,'*');
    end
end

[~,bF]=histc(frame,linspace(1,9000,5));
[~,bS]=histc(sess,cat(2,1,quantile(sess,10),Inf));
bdT = cat(2,0,quantile(dnT(k),30),Inf);
[~,bdnT]=histc(dnT,bdT);


% construct t2s and hitmiss for nF
hitmiss1 = [];
t2s1 = [];
for iM=1:nM
    for iS=1:nS
        hitmiss1 = cat(1,hitmiss1,repmat(hitmiss(iM,iS),[nframe(iM,iS),1]));
        t2s1 = cat(1,t2s1,repmat(t2s(iM,iS),[nframe(iM,iS),1]));
    end
end
k = nV>0.1 & ~edge & ~isnan(znidphi1) & dnT<15 & t2s1 > -40;
[p,table,stats] = anovan(znidphi1(k),{bait(k),conc(k),hitmiss1(k),trial(k)},'varnames',{'bait','conc','hitmiss','trial'},'model','interaction');


k = nV>0.1 & ~edge & ~isnan(znidphi1) & t2s1 > -40;
H = histcn([dnT(k),hitmiss1(k)],bdT,1:2,'AccumData',znidphi1(k),'fun',@nanmedian);
dH = histcn([dnT(k),hitmiss1(k)],bdT,1:2,'AccumData',znidphi1(k),'fun',@nanstderr);
N = histcn([dnT(k),hitmiss1(k)],bdT,1:2);
clf;
lineProps.col = {'b'};
mseb(bdT(1:end-1),H(:,1)',dH(:,1)',lineProps);
lineProps.col = {'r'};
mseb(bdT(1:end-1),H(:,2)',dH(:,2)',lineProps);
line([0 90],[0 0],'color','k');
legend('found','missed');

% Conclusion: KP sessions mice run ~5cm/s faster.  Could be because they
% are dumped in via bucket instead of placed on by hand

%
%
%
timeInRing


hitmiss2 = permute(repmat(hitmiss,[1,1,31]),[3,1,2]);
t2sx = permute(repmat(t2s,[1,1,31]),[3,1,2]);

k = ~isnan(R) & t2sx > -40;
[p,table,stats]=anovan(R(k),{radii0(k),mouse0(k),conc0(k),bait0(k),ALorKP0(k),trial0(k),hitmiss2(k)},'varnames',{'radius','mouse','concentration','bait','ALorKP','trial','spotfound'},'model','interaction');


clf;
lineProps.col = {'b'};
kfound = t2s>0;
knot = t2s ==-20;
%ALorKP0(isnan(ALorKP0))=0;
k1 = squeeze(ALorKP0(1,:,:));% & kfound;
mseb(radii,nanmean(R(:,k1==0),2),nanstderr(R(:,k1==0),[],2)',lineProps);
lineProps.col = {'r'};
mseb(radii,nanmean(R(:,k1==1),2),nanstderr(R(:,k1==1),[],2)',lineProps);
lineProps.col = {[0.5 0.5 0.5]};
mseb(radii,nanmean(randR(:,:),2),nanstderr(randR(:,:),[],2)',lineProps);
line([0 30],[0 0],'color','k');
legend('AL','KP','rand');

clf;
plot(radii,log10(nanmean(R(:,k1==0),2)),'r','linewidth',2);
hold on;
plot(radii,log10(nanmean(R(:,k1==1),2)),'b','linewidth',2);
plot(radii,log10(nanmean(randR(:,:),2)),'color',[0.5 0.5 0.5],'linewidth',2);
%line([0 30],[0 0],'color','k');
legend('AL','KP','rand');

clf;
lineProps.col = {'b'};
kfound = t2s>0;
knot = t2s ==-20;
%ALorKP0(isnan(ALorKP0))=0;
%k1 = squeeze(ALorKP0(1,:,:));% & kfound;
mseb(radii,nanmean(R(:,kfound),2),nanstderr(R(:,k1==0),[],2)',lineProps);
lineProps.col = {'r'};
mseb(radii,nanmean(R(:,knot),2),nanstderr(R(:,k1==1),[],2)',lineProps);
lineProps.col = {[0.5 0.5 0.5]};
mseb(radii,nanmean(randR(:,:),2),nanstderr(randR(:,:),[],2)',lineProps);
line([0 30],[0 0],'color','k');
legend('found','missed','rand');


clf;
plot(radii,log10(nanmean(R(:,kfound==1),2)),'r','linewidth',2);
hold on;
plot(radii,log10(nanmean(R(:,knot==1),2)),'b','linewidth',2);
plot(radii,log10(nanmean(randR(:,:),2)),'color',[0.5 0.5 0.5],'linewidth',2);
%line([0 30],[0 0],'color','k');
legend('found','missed','rand');


k = ~isnan(mouse0);
kfound = t2s>0;
knot = t2s ==-20;
clear p h stats
clear p2 h2 stats2

for iR=2:size(R,1)
    disp(iR);
    [p(iR),h(iR)] = ranksum(R(iR,k & kfound),R(iR,k & knot));
end

clf;
lineProps.col = {'b'};
mseb(radii,nanmedian(R(:,kfound),2),nanstderr(R(:,kfound),[],2)',lineProps);
lineProps.col = {'r'};
mseb(radii,nanmedian(R(:,knot),2),nanstderr(R(:,knot),[],2)',lineProps);
lineProps.col = {[0.5 0.5 0.5]};
mseb(radii,nanmedian(randR(:,:),2),nanstderr(randR(:,:),[],2)',lineProps);
line([0 30],[0 0],'color','k');

expfalsePos = size(R,1);   
for iR=1:size(R,1)
    if p(iR)<0.05/expfalsePos
       text(radii(iR),nanmean(R(iR,kfound),2)+1E-5,'*');
    end
end

legend('spot found','spot not found','random')
xlabel('distance from spot (cm)','fontsize',18);
ylabel('time / area (s/cm^2)','fontsize',18);


% line props arguments
blue.col = {'b'};
black.col = {'k'};
red.col = {'r'};
green.col = {'g'};
cyan.col = {'c'};

[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),2),Inf));
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
sess1 = permute(repmat(bS,[1,1,31]),[3,1,2]);
sess1(sess1==0)=nan;
spotfound1 = permute(repmat(spotfound,[1,1,31]),[3,1,2]);
bait1 = permute(repmat(bait0,[1,1,31]),[3,1,2]);
trial1 = permute(repmat(trial0,[1,1,31]),[3,1,2]);
mouse1 = permute(repmat(mouse0,[1,1,31]),[3,1,2]);
initd1 = permute(repmat(initD,[1,1,31]),[3,1,2]);


% test for variability due to mouse
k = ~isnan(Rrat) & Rrat>0;
[p,table,stats] = anovan(log10(Rrat(k)),{mouse1(k)},'varnames',{'mouse'},'model','full');

% remove variability due to mouse
zRrat = nan(size(R)); for iM=1:nM; k = mouse1==iM & Rrat~=0; zRrat(k) = nanzscore(log10(Rrat(k))); end;
zR = nan(size(R)); for iM=1:nM; k = mouse1==iM & R~=0; zR(k) = nanzscore(log10(R(k))); end;
[p,table,stats] = anovan(zRrat(:),{mouse1(:)},'varnames',{'mouse'},'model','full');

% test effect of session and trial
[p,table,stats] = anovan(zRrat(:),{trial1(:),sess1(:)},'varnames',{'trial','sess1'},'continuous',2,'model','full');

% for some reason initial distance for even trials is lower (?)
[p,table,stats]=anova1(initD(:),trial0(:));

% for some reason even / odd trials group together (?)
clf;
mseb(radii,(nanmean((zRrat(:,trial0==1)),2))',nanstderr(zRrat(:,trial0==1),[],2)',blue);
mseb(radii,(nanmean((zRrat(:,trial0==2)),2))',nanstderr(zRrat(:,trial0==2),[],2)',black);
mseb(radii,(nanmean((zRrat(:,trial0==3)),2))',nanstderr(zRrat(:,trial0==3),[],2)',green);
mseb(radii,(nanmean((zRrat(:,trial0==4)),2))',nanstderr(zRrat(:,trial0==4),[],2)',red);
mseb(radii,(nanmean((zRrat(:,trial0==5)),2))',nanstderr(zRrat(:,trial0==5),[],2)',cyan);
legend('1','2','3','4','5');

% test for effect of radius, bait, spot found
[p,table,stats] = anovan(zRrat(:),{radii1(:),bait1(:),spotfound1(:)},'varnames',{'radii','bait','found'},'model','full');


% limit to radius <= 15cm
k = radii1 <= 15;

[p,table,stats] = anovan(zRrat(k),{bait1(k),spotfound1(k)},'varnames',{'bait','found'},'model','full');
clf;
mseb(radii,(nanmean((zRrat(:,bait0==0 & spotfound==0)),2))',nanstderr(zRrat(:,bait0==0 & spotfound==0),[],2)',blue);
mseb(radii,(nanmean((zRrat(:,bait0==0 & spotfound==1)),2))',nanstderr(zRrat(:,bait0==0 & spotfound==1),[],2)',black);
mseb(radii,(nanmean((zRrat(:,bait0==1 & spotfound==0)),2))',nanstderr(zRrat(:,bait0==1 & spotfound==0),[],2)',green);
mseb(radii,(nanmean((zRrat(:,bait0==1 & spotfound==1)),2))',nanstderr(zRrat(:,bait0==1 & spotfound==1),[],2)',red);
legend('not-not','not-found','bait-not','bait-found');

% test for effect of bait and prev trial...some minor interactions
[p,table,stats] = anovan(zRrat(k),{prevT1(k),bait1(k),spotfound1(k)},'varnames',{'prevT','bait','found'},'model','full');
boxplot(zRrat(k),{prevT1(k)},'notch','on');

% mice explore closer to the spot more frequently following unbaited trials

clf;
mseb(radii,(nanmean((zRrat(:,prevT0==0)),2))',nanstderr(zRrat(:,prevT0==0),[],2)',blue);
mseb(radii,(nanmean((zRrat(:,prevT0==1)),2))',nanstderr(zRrat(:,prevT0==1),[],2)',black);
legend('prev trial unbaited','prev trial baited');

clf;
mseb(radii,(nanmean((zRrat(:,prevT0==0 & bait0==0)),2))',nanstderr(zRrat(:,prevT0==0 & bait0==0),[],2)',blue);
mseb(radii,(nanmean((zRrat(:,prevT0==0 & bait0==1)),2))',nanstderr(zRrat(:,prevT0==0 & bait0==1),[],2)',black);
mseb(radii,(nanmean((zRrat(:,prevT0==1 & bait0==0)),2))',nanstderr(zRrat(:,prevT0==1 & bait0==0),[],2)',green);
mseb(radii,(nanmean((zRrat(:,prevT0==1 & bait0==1)),2))',nanstderr(zRrat(:,prevT0==1 & bait0==1),[],2)',red);
legend('unbaited-unbaited','unbaited-baited','baited-unbaited','baited-baited');

