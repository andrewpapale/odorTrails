% 2018-11-19 AndyP
% spot stats

%% Basic Stats

load('2019-07-08-FullDataset.mat','mint2s','edge','isstopped','edgespot','t2s0','spotfound0trunc','removeint2s','nV','V','spotfound','t2s');

mint2s = 1;

% number of frames for most analyses
totframes = sum(~edge & ~isstopped & ~edgespot & (t2s0==1 | spotfound0trunc==1) & ~removeint2s & nV<100 & V<100);
fprintf('total frames for analysis = %0.2f \n',totframes);

% % of frames for most analyses
fracframes = sum(~edge & ~isstopped & ~edgespot & (t2s0==1 | spotfound0trunc==1) & ~removeint2s & nV<100 & V<100)./length(edge);
fprintf('% frames for analysis = %0.2f \n',fracframes*100);

% number of frames for Figure2C
framesexplo = sum(~edge & ~edgespot & ~removeint2s & nV<100 & V<100);
fprintf('number of frames for Figure2C = %0.2f \n',framesexplo);

nS = sum((spotfound(:)==1 & t2s(:)>mint2s) | (t2s(:)>-40 & spotfound(:)==0)); % number of sessions
nSucc = sum(spotfound(:)==1 & t2s(:)>mint2s); % number of successful trials
fracSucc = nSucc./nS; % fraction of successful trials
fprintf('total number of trials = %0.2f \n',nS);
fprintf('number of successful trials = %0.2f \n',nSucc);
fprintf('% of successful trials = %0.2f \n',fracSucc*100);

load('2019-08-12-shuffleSpots.mat','sht2s');
spotfound = nan(size(sht2s));
spotfound(sht2s > mint2s)=1;
spotfound(sht2s==-20)=0;
nS = nansum(spotfound(:)==1 | spotfound(:)==0);
nSucc = nansum(spotfound(:)==1);
fracSucc = nSucc./nS*100;

%% Logistic Regression on Successful Trials
% statsFigure2A
% 2018-09-14 AndyP
% logistic regression on % correct

load('2019-07-08-FullDataset.mat','Dnid','idphi','tstopped','isstopped','edge','edgespot','t2s0','spotfound0trunc','removeint2s','nV','V','t2s','mint2s','spotfound','sess0','pExplored','ALorKP0','bait0','conc0','trial0','mouse0','initD','initO','meanV','meannV','fractoward','ddnV','mouse','sess','dnT');
k = dnT < 30 & ~isstopped & ~edge & (t2s0==1 | spotfound1==0) & nV < 100 & V<100 & ~removeint2s & ~edgespot;
H1 = histcn([mouse(k),sess(k)],1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N1 = histcn([mouse(k),sess(k)],1:5,1:197);
H1(N1==0)=nan;
H1(:,end+1)=nan;

mint2s = 1;
k = (t2s>mint2s) | (spotfound==0 & t2s~=-40);

[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf)); % bin session
bS = reshape(bS,[size(sess0,1),size(sess0,2)]);
bS(bS==0)=nan;
pE = squeeze(pExplored(12,:,:));

X1 = [ALorKP0(k),bait0(k),conc0(k),trial0(k),bS(k),mouse0(k),initD(k),initO(k),meanV(k),meannV(k),Dnid(k),H1(k),pE(k),ddnV(k),tstopped(k)];
Xvars = {'','ALorKP','bait','conc','trial','sess','mouse','initD','initO','meanV','meannV','disprat','meandphi','pE','ddnV','tstopped'}';
Y = spotfound(k);
options.MaxIter = 10000;
[b0,~,stats] = glmfit(X1,Y,'binomial','link','logit','options',options);
table(Xvars,stats.p, b0, stats.p<0.05/15)

%% % correct data vs shuffles

load('2019-07-08-FullDataset.mat','t2s','spotfound','mouse0');

mint2s = 1;

k = t2s~=-40 & (t2s==-20 | t2s>mint2s & spotfound==1); 
DATA = histcn(mouse0(k),1:5,'AccumData',spotfound(k),'fun',@nanmean);

load('2019-08-12-shuffleSpots.mat','sht2s');
ish = permute(repmat((1:100)',[1,5,197]),[2,3,1]);
spotfound = nan(size(sht2s));
spotfound(sht2s>1)=1;
spotfound(sht2s==-20)=0;

SHUFF = histcn(ish(:),1:100,'AccumData',spotfound(:),'fun',@nanmean);

[h,p,ci,stats] = ttest2(DATA,SHUFF,'Vartype','Unequal','Tail','Right');
fprintf('t-test mouse vs. shuffles p=%0.2e,t-stat=%0.2f,df=%0.2f \n',p,stats.tstat,stats.df);
    


%% Time to Spot
load('2019-07-08-FullDataset.mat','t2s','mint2s','initD');
k = t2s>mint2s; % median time to spot for successful trials
DATA = t2s(k)./initD(k);
mt2s = nanmean(DATA);
fprintf('time to spot / init distance (s/cm) MICE = %0.2f \n',mt2s);

load('2019-08-12-shuffleSpots.mat','sht2s','sinitD');
k = sht2s>1;
SHUFF = sht2s(k)./(sinitD(k)./11.2);
smt2s = nanmean(SHUFF);
fprintf('time to spot / init distance (s/cm) SHUFF = %0.2f \n',smt2s);


[h,p,ci,stats]=ttest2(DATA,SHUFF,'Vartype','unequal','Tail','Left');
fprintf('t-test mouse vs. CRW p=%0.2e,t-stat=%0.2f,df=%0.2f \n',p,stats.tstat,stats.df);


%% Displacement from initial Distance
load('2019-07-08-FullDataset.mat','Dnid','mint2s','initD');
k = t2s>mint2s; % median time to spot for successful trials
DATA = Dnid(k);
md = nanmean(DATA);
fprintf('nose displacement / init distance (cm/cm) MICE = %0.2f \n',md);

load('2019-08-12-shuffleSpots.mat','sht2s','sinitD','sDnid');
k = sht2s>1;
SHUFF = sDnid(k).*11.2;
smd = nanmean(SHUFF);
fprintf('nose displacement / init distance (cm/cm) SHUFF = %0.2f \n',smd);


[h,p,ci,stats]=ttest2(DATA,SHUFF,'Vartype','unequal','Tail','Left');
fprintf('t-test mouse vs. CRW p=%0.2e,t-stat=%0.2f,df=%0.2f \n',p,stats.tstat,stats.df);

%% Partitioning Table Analysis

load('2019-07-08-FullDataset.mat','mint2s','parts','pExplored','t2s','sess0','spotfound');
maxX = 1280;
maxY = 1020;
nP = length(parts);
bcmx = nan(size(parts));
bcmy = nan(size(parts));
for iN=1:nP
    nbinsx = linspace(44.8,maxX-44.8,ceil(sqrt(parts(iN))));
    nbinsy = linspace(44.8,maxY-44.8,ceil(sqrt(parts(iN))));
    bcmx(iN) = nanmedian(diff(nbinsx))./11.2;
    bcmy(iN) = nanmedian(diff(nbinsy))./11.2;
end
bcmxy = sqrt(bcmx.^2+bcmy.^2);

% percent of table explored
pE = squeeze(pExplored(12,:,:));
nanmean(pE(t2s>mint2s))
nanstd(pE(t2s>mint2s))

db = repmat(bcmxy',[1,size(pExplored,2),size(pExplored,3)]);
df = permute(repmat(spotfound,[1,1,size(bcmxy,2)]),[3,1,2]);

[p,tab,stats]=anovan(pExplored(:),{db(:),df(:)},'model','full');

[~,bS]=histc(sess0(:),cat(2,1,quantile(sess0(:),5),Inf)); % bin session
bS(bS==0)=nan;
pE = squeeze(pExplored(12,:,:));
[p,tab,stats] = anova1(pE(:),bS(:));

%% Occupancy Analysis
load('2019-07-08-FullDataset.mat','R','mint2s','t2s','area','radii');
% we want different variables
mint2s = 1;
k = t2s>mint2s;
pRs = R(:,k);
N = nansum(pRs,2);
pRs(N<30,:)=nan;
N0 = nansum(pRs,1);
N0(N0<30)=nan;
mpRsD = pRs./(repmat(N0,[size(pRs,1),1]).*area);
k = t2s==-20;
pRn = R(:,k);
N = nansum(pRn,2);
pRn(N<30,:)=nan;
N0 = nansum(pRn,1);
N0(N0<30)=nan;
mpRC = pRn./(repmat(N0,[size(pRn,1),1]).*area);

nR = length(radii);
p = nan(nR,1);
kstat = nan(nR,1);
for iR=1:nR
    try
        [~,p(iR),kstat0]=kstest2(mpRsD(iR,:),mpRC(iR,:),'Tail','smaller');
        kstat(iR)=kstat0;
    catch
    end
end

hbon = p < 0.05./nR;

table(radii',hbon,p,kstat)
%% Velocity


load('2019-07-08-FullDataset.mat','spotfound','mouse0','sess0','nM','nS','removeint2s','edge','isstopped','t2s0','spotfound0trunc','edgespot','V','nV','dnT','spotfound1','mouse','sess');

bins = linspace(1.5,100,30);
distmat = repmat(bins(1:end-1)',[1,5,197]);
distmat = distmat(:,:);

% nose
k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;
meannV = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197,'AccumData',nV(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197);
meannV(N==0)=nan;

mH = squeeze(nanmean(meannV(:,:),2));
dH = squeeze(nanstderr(meannV(:,:),[],2));

[p,tab,stats]=anova1(meannV(:),distmat(:));

% body
k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;
meanV = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197,'AccumData',V(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197);
meanV(N==0)=nan;

mH = squeeze(nanmean(meanV(:,:),2));

dH = squeeze(nanstderr(meanV(:,:),[],2));
[p,tab,stats]=anova1(meanV(:),distmat(:));

%%  orientation

load('2019-07-08-FullDataset.mat','spotfound0trunc','orient','dnT','t2s0','edge','isstopped','edgespot','removeint2s','V','nV');

% uses circular statistics toolbox...
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
k = ~edge & ~isstopped & t2s0==1 & ~edgespot & ~removeint2s & V < 100 & nV < 100;
[~,bins,bdnT]=histcounts(dnT,linspace(1.5,100,50));
bdnT(bdnT==0)=nan;
k0 = k & ~isnan(bdnT) & ~isnan(orient);
pval=circ_wwtest(orient(k0)*pi/180,bdnT(k0));
%
pval1 = nan(49,1);
for iD=1:49
    k0 = k & bdnT==iD;
    [pval1(iD)]=circ_medtest(orient(k0)*pi/180,0);
end

H = pval1<0.05/49;
table(bins(1:49)',pval1,H)

%% casting curvature


bins = linspace(1.5,100,30);


log10dphi = logVar(abs(nidphi)./nV);

k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;

meanC = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
meanC(N==0)=nan;

group1 = size(meanC(:,:),2);

distmat = repmat(bins(1:end-1)',[1,group1]);

[p,tab,stats]=anova1(meanC(:),distmat(:));

meanC = histcn([dnT(k),spotfound1(k),mouse(k),sess(k)],bins,0:1,1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([dnT(k),spotfound1(k),mouse(k),sess(k)],bins,0:1,1:5,1:197);
meanC(N==0)=nan;

mHf = squeeze(meanC(:,2,:));
mHn = squeeze(meanC(:,1,:));

group1 = ones(size(mHf));
group2 = ones(size(mHn))+1;

[p,tab,stats]=anovan(cat(1,mHf(:),mHn(:)),{cat(1,distmat(:),distmat(:)),cat(1,group1(:),group2(:))},'model','full');

% curvature vs session
k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;
mH = histcn([mouse(k),sess(k)],1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([mouse(k),sess(k)],1:5,1:197);
mH(N==0)=nan;
mH = mH(:,:);
bS = repmat((1:197),[5,1]);
bS = bS(:);
[~,~,bS]=histcounts(bS,6);

[p,tab,stats]=anova1(mH(:),bS(:));

% curvature vs velocity or distance GLM
k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;

meanC = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
meanC(N==0)=nan;

meanV = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',log10(nV(k)),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
meanV(N==0)=nan;
meanV(isinf(meanV))=nan;

X = [distmat(:),meanV(:)];
Y = meanC(:);
stepwiseglm(X,Y,'constant');

% nose velocity vs session
k = ~edge & ~isstopped & (t2s0==1 | spotfound1==0) & ~edgespot & ~removeint2s & nV < 100 & V < 100;
mH = histcn([mouse(k),sess(k)],1:5,1:197,'AccumData',nV(k),'fun',@nanmean);
N = histcn([mouse(k),sess(k)],1:5,1:197);
mH(N==0)=nan;
mH = mH(:,:);
bS = repmat((1:197),[5,1]);
bS = bS(:);
[~,~,bS]=histcounts(bS,6);

[p,tab,stats]=anova1(mH(:),bS(:));
%%
load('full_model_30_seconds_1.5cm_02012019_vars.mat','cx','cy','sess','spotfound','x','y','nx','ny','dnT');

morient = nan(size(x));
nS = length(spotfound);
for iS=1:nS
    k = sess==iS;
    A = atan2d(y(k)-cy(iS),x(k)-cx(iS));
    B = atan2d(y(k)-ny(k),x(k)-nx(k));
    morient(k) = ((A)-(B));
    disp(iS);
end




figure(1); clf;
H = histogram2(dnT(:),(morient(:)),linspace(0,100,50),linspace(-180,180,45),'facecolor','flat');

% k = ~edge & ~isstopped & spotfound0trunc==1;
% F2 = figure(2); clf;
% Hnot = histogram2(dnT(k),abs(orient(k)),linspace(0,120,120),linspace(0,180,45),'facecolor','flat');

% H1 = Hnot.Values;

H0 = H.Values;
Hx = H.XBinEdges;
Hy = H.YBinEdges;


nH = size(H0,1);
for iH=1:nH
    H0(iH,:)=H0(iH,:)./nansum(H0(iH,:));
    %H1(iH,:)=H1(iH,:)./nansum(H1(iH,:));
end
H0(H.Values==0)=nan;

F = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial'); 
clf;
%pcolor(Hx(1:size(H0,1)),Hy(1:size(H0,2)),smooth2a(H0'-H1',2,2)); shading flat;
pcolor(Hx(1:size(H0,1)),Hy(1:size(H0,2)),smooth2a(H0',2,2)); shading flat;




[~,bins,bdnT]=histcounts(dnT,linspace(1.5,100,50));
bdnT(bdnT==0)=nan;
k1 = ~isnan(bdnT) & ~isnan(morient);
pval=circ_wwtest(morient(k1)*pi/180,bdnT(k1));


pval1 = nan(49,1);
for iD=1:49
    k0 = k1 & bdnT==iD;
    [pval1(iD)]=circ_medtest(morient(k0)*pi/180,circ_ang2rad(0));
end

H = pval1<0.05;
table(bins(1:49)',pval1,H)

%%
pval1 = nan(49,1);
for iD=1:49
    k0 = bdnT==iD & ~isnan(morient);
    %[pval1(iD)]=circ_medtest(morient(k0)*pi/180,circ_ang2rad(0));
    [h(iD) mu(iD)] = circ_mtest(morient(k0)*pi/180, circ_ang2rad(0), 0.05/49);
end

%H = pval1<0.05/49;
%table(bins(1:49)',pval1,H)
%%
% k = t2s0==0;
% F2 = figure(2); clf;
% Hnot = histogram2(dnTC(k),abs(rOC(k)),linspace(0,120,120),linspace(0,180,180),'facecolor','flat');

%H1 = Hnot.Values;

H0 = H.Values;
Hx = H.XBinEdges;
Hy = H.YBinEdges;



nH = size(H0,1);
for iH=1:nH
    HC(iH,:)=H0(iH,:)./nansum(H0(iH,:));
    %     H1(iH,:)=H1(iH,:)./nansum(H1(iH,:));
end
HC(H.Values==0)=nan;

F = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial'); clf;
pcolor(Hx(1:size(H0,1)),Hy(1:size(H0,2)),smooth2a(HC',2,2)); shading flat;
C = colorbar;
set(C,'fontsize',36,'FontName','Arial');
set(C,'Limits',[0 0.05]);
set(C,'Ticks',[0 0.05]);
caxis([0 0.05]);
set(gca,'XLim',[0 100]);
set(gca,'fontsize',36);
set(gca,'YTick',[0 90 180],'YLim',[0 180]);
box off;
xlabel('Distance from Spot (cm)','fontsize',24,'FontName','Arial');
ylabel(sprintf('angle relative to spot (%c)', char(176)),'fontsize',24,'FontName','Arial');

% radii = linspace(0,120,120);
% p=nan(size(HC,1),1);
% k=nan(size(HC,1),1);
% for iR=1:size(HC,1)
%     k1 = ~isnan(HC(iR,:));
%     k2 = ~isnan(HD(iR,:));
%     A = HC(iR,k1);
%     B = HD(iR,k2);
%     if ~isempty(A) && ~isempty(B)
%         [~,p(iR),k(iR)]=kstest2(A(:),B(:));
%     end
% end
% 
% h = p < 0.05/size(HC,1);
% 
% table(radii(1:size(p,1))',h,p,k)





%% Casting Distance Measure


% polar_phi vs session
minB = 1;
nB = 10;
bins = linspace(1.5,100,nB);
k = ~edge & ~isstopped & t2s0==1 & ~edgespot & ~removeint2s & nV < 100 & V < 100;
mH = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',abs(polar_phi(k))./(360*dnT(k)),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
mH(N<minB)=nan;
mH = mH./N;
mH = mH(:,:);
mH = mH(1,:);
bS = repmat((1:196),[5,1]);
bS = bS(:);
[~,~,bS]=histcounts(bS,6);

[p,tab,stats]=anova1(mH,bS(:));

% curvature vs trial order
minB = 1;
bins = linspace(0,100,25);
k = ~edge & ~isstopped & t2s0==1 & ~edgespot & ~removeint2s & nV < 100 & V < 100;
mH = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
mH(N<minB)=nan;
mH = mH(:,:);
mH = mH(1,:);
group = trial0(:,1:196);
group = group(:);
group(group==5)=nan;
group2 = bait0(:,1:196);
[p,tab,stats] = anova1(mH,group);
[p,tab,stats] = anovan(mH,{group,group2(:)},'model','full');

% nose deflection vs trial order
minB = 1;
bins = linspace(0,100,25);
k = ~edge & ~isstopped & t2s0==1 & ~edgespot & ~removeint2s & nV < 100 & V < 100;
mH = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197,'AccumData',log10(abs(Dline(k))),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:5,1:197);
mH(N<minB)=nan;
mH = mH(:,:);
mH = mH(1,:);
group = trial0(:,1:196);
group = group(:);
group(group==5)=nan;
[p,tab,stats] = anova1(mH,group);

%% Model stats

% logistic regression on % correct for parameter screen

load('2019-08-26-50k param sweep.mat');

X1 = cat(2,Vmax,kbin,kint,k,noiselevel,filt);
Y = successes;
options.MaxIter = 10000;
[b0,~,stats] = glmfit(X1,Y,'binomial','link','logit','options',options);
Xvars = {'','Vmax','kbin','kint','k','noiselevel','filt'};
table(Xvars',stats.p, b0, stats.p<0.05/6)

% Figure3F Reduced Models
load('par_screen_success_rates_1000_reps_1.5cm_02032019.mat');


for iD=1:length(success_data)
    switch success_data(iD).model
        case 'VSC'
            group0 = 1;
        case 'SC'
            group0 = 2;
        case 'VC'
            group0 = 3;
        case 'VS'
            group0 = 4;
        case 'V'
            group0 = 5;
        case 'S'
            group0 = 6;
        case 'C'
            group0 = 7;
        case 'null'
            group0 = 8;
    end
    group(iD) = group0;
    
    success(iD) = success_data(iD).successrate;
end

% [p,tab,stats]=anova1(success,group);
% c = multcompare(stats,'ctype','bonferroni');


load('CRW_all_sims_1.5cm_02012019_vars.mat', 'spotfound','type');
CRWsucc = histcn(type,1:8,'AccumData',spotfound,'fun',@nanmean);

succ1 = cat(2,success,CRWsucc');
group1 = cat(2,group,repmat(9,[1,8]));
[p,tab,stats]=anova1(succ1,group1);
c = multcompare(stats,'ctype','bonferroni');

%% Model Occupancy

%% Occupancy Analysis

radius = linspace(1.5,100,100);
nR = length(radius);
area = nan(nR,1);
for iR=1:nR
    if iR==1
        area(iR)=pi*radius(iR).^2;
    else
        area(iR)=pi*(radius(iR).^2-radius(iR-1).^2);
    end
end
radii = radius;


load('full_model_30_seconds_1.5cm_kc_1_sigfilt_4_30cmthreshold_09132019_vars.mat','Rm','t2s');

mint2s = 1;
minB = 1;
%area = 1;
k = t2s>mint2s;
pRs = Rm(:,k);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
N0 = nansum(pRs,1);
N0(N0<minB)=nan;
mpRsCSM = pRs./(repmat(N0,[size(pRs,1),1]).*area);



load('CRW_all_sims_1.5cm_08152019_vars.mat','Rc','type');

% 
%area = 1;
pRs = Rc(:,:);
N = nansum(pRs,2);
pRs(N<minB,:)=nan;
N0 = nansum(pRs,1);
N0(N0<minB)=nan;
mpRsCRW = pRs./(repmat(N0,[size(pRs,1),1]).*area);

nR = length(radii);
p = nan(nR,1);
kstat = nan(nR,1);
for iR=1:nR
    try
        [~,p(iR),kstat0]=kstest2(mpRsCSM(iR,:),mpRsCRW(iR,:),'Tail','smaller');
        kstat(iR)=kstat0;
    catch
    end
end

hbon = p < 0.05./nR;

table(radii',hbon,p,kstat)

