% plot random sessions

nM = max(mouse);
nS = max(sess);

i=-1;
while i<1
    rm = randi(nM,1);
    rs = randi(nS,1);
    
    k = mouse==rm & sess==rs;
    
    if sum(k)>0
        plot(x(k),y(k),'k.');
        set(gca,'XLim',[1 1280],'YLim',[1 1040]);
        pause;  
    end
end



for iM=1:nM
    rs = randperm(nS,nS);
    for iS=1:nS
        
        k = mouse==iM & sess==rs(iS);
        
        if sum(k)>0
            %plot(x(k),y(k),'k.');
            scatter(x(k),y(k),20,V(k),'filled');
            set(gca,'XLim',[1 1280],'YLim',[1 1040]);
            caxis([0 100]);
            title(iM);
            pause;
        end
    end
end
    