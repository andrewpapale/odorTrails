bins = linspace(1.5,100,30);
log10dphi = logVar(abs(ndphi)./nV);

k = ~edge & t2s2==1 & nV < 100 & bV < 100;


meanC = histcn([dnT(k),cond0(k),mouse(k),sess(k)],bins,2:3,1:5,1:197,'AccumData',log10dphi(k),'fun',@nanmean);
N = histcn([dnT(k),cond0(k),mouse(k),sess(k)],bins,2:3,1:5,1:197);
meanC(N==0)=nan;

mHf = squeeze(nanmean(meanC(:,2,:),3));
dHf = squeeze(nanstderr(meanC(:,2,:),[],3));
mHn = squeeze(nanmean(meanC(:,1,:),3));
dHn = squeeze(nanstderr(meanC(:,1,:),[],3));

F = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial');
clf;
lineProps.col = {'r'};
mseb(bins(1:length(mHf)),mHf,dHf',lineProps);
lineProps.col = {'b'};
mseb(bins(1:length(mHn)),mHn,dHn',lineProps);
set(gca,'fontsize',36);
xlabel('Distance from Spot (cm)','fontsize',24);
ylabel('log_{10}(|dphi|/nV) [cm^{-1}]','fontsize',24);
L = line([1.5 1.5],[-1.4 -0.6]);
set(L,'LineStyle','--','color','k','LineWidth',2);
L = legend('Spot Found','Not Found');
set(L,'fontsize',24);
%set(gca,'YLim',[-1.4 -0.6]);
