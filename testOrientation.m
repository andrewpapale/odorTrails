

k = t2s0==1 & bV < 100 & nV < 100 & t2s2==1 & cond0==0;
figure(1); clf;
H = histogram2(dnT(k),abs(orient(k)),linspace(0,100,50),linspace(0,180,45),'facecolor','flat');

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
C = colorbar;
set(C,'fontsize',21);
set(C,'Limits',[0 0.05]);
set(C,'Ticks',[0 0.05],'TickLabels',{'0','5'});
caxis([0 0.05]);
set(gca,'XLim',[0 100]);
set(gca,'fontsize',24);
set(gca,'YTick',[0 90 180],'YLim',[0 180]);
box off;
xlabel('Distance from Spot (cm)','fontsize',24);
ylabel(sprintf('angle relative to spot (%c)', char(176)),'fontsize',24);
