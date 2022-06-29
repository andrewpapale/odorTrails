iC = 1;
xbins = linspace(0,600,60); 
ybins = linspace(0,500,40);
[X,Y]=meshgrid(xbins(1:end-1),ybins(1:end-1));
X = X';
Y = Y';
figure(1); clf;
for iP=1:4
    for iW=0:1
        
        k = ~edge & ~ktoss & confn > 0.3 & t2s2 & sP0==iP & wind0==iW;
        H = histcn([x(k),y(k)],xbins,ybins,'AccumData',tro(k)*pi/180,'fun',@circ_mean);
        dH = histcn([x(k),y(k)],xbins,ybins,'AccumData',tro(k)*pi/180,'fun',@circ_std);
        dH(dH==0)=nan;
        N = histcn([x(k),y(k)],xbins,ybins);
        H(N==0)=nan;
        
        subplot(4,2,iC);
        k0 = abs(sin(H(:))./dH(:)) < pi & abs(cos(H(:))./dH(:)) < pi;
        quiver(X(k0)./5.174,Y(k0)./5.174,cos(H(k0))./dH(k0),sin(H(k0))./dH(k0),2,'LineWidth',1.5,'MaxHeadSize',20,'AlignVertexCenters','on');
        hold on;
        viscircles([cx0(sP==iP & ~removeS)./5.174,cy0(sP==iP & ~removeS)./5.174],repmat(1.5,[sum(sP(~removeS)==iP),1]));
        
        iC = iC+1;
        
    end
end