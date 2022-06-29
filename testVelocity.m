bins = linspace(1.5,100,30);

k = ~edge & t2s2==1 & nV < 100 & bV < 100 & cond0==1;
meannV = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197,'AccumData',nV(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197);
meannV(N==0)=nan;

mHf = squeeze(nanmean(meannV(:,:),2));
dHf = squeeze(nanstderr(meannV(:,:),[],2));

k = ~edge & cond0==0 & t2s2==1 & nV < 100 & bV < 100;
meannV = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197,'AccumData',nV(k),'fun',@nanmean);
N = histcn([dnT(k),mouse(k),sess(k)],bins,1:6,1:197);
meannV(N==0)=nan;

mHn = squeeze(nanmean(meannV(:,:),2));
dHn = squeeze(nanstderr(meannV(:,:),[],2));


F = figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial');
clf;
lineProps.col = {'r'};
mseb(bins(1:length(mHf)),mHf,dHf',lineProps);
lineProps.col = {'b'};
mseb(bins(1:length(mHn)),mHn,dHn',lineProps);
set(gca,'fontsize',36);
xlabel('Distance from Spot (cm)','fontsize',24);
ylabel('Mean Nose Velocity (cm/s)','fontsize',24);
L = line([1.5 1.5],[0 100]);
set(L,'LineStyle','--','color','k','LineWidth',2);
L = legend('Spot Found','Not Found');
set(L,'fontsize',24);
%set(gca,'YLim',[12 30]);