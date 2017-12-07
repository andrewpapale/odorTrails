% plot spots
nM = max(mouse);
nS = max(sess);
F = Figure(1); clf;
for iM=1:nM
    for iS=1:nS
        k = mouseT==iM & sessT==iS;
        
        xT = xT1(k);
        yT = yT1(k);
        
        color = [0 0 0];
        while sum(color)<0.2
            color = rand(3,1);
        end
            
        
        plot(xT,yT,'.','color',color);
        hold on;
        
        xT0 = nanmedian(xT);
        yT0 = nanmedian(yT);
        
        if xT0<44.8 | xT0>1280-44.8 | yT0<44.8 | yT0>1024-44.8 %#ok<OR2>
            plot(xT,yT,'.','color',[0 0 0]);
        end
    end
end
    
    set(gca,'XLim',[1 1280],'YLim',[1 1024]);
    axis off;
    line([44.8 44.8],[1 1024],'color','r');
    line([1280-44.8 1280-44.8],[1 1024],'color','r');
    line([1 1280],[44.8 44.8],'color','r');
    line([1 1280],[1024-44.8 1024-44.8],'color','r');
    
F = figure(2);
getTime2Spot;
k = t2s>-40 & ~isnan(mouse0);
scatterhist(xT(k),yT(k),'group',ALorKP0(k),'Kernel','on');

