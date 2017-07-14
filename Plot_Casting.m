multipleSess = true;


if multipleSess
    nM = length(unique(mouse));
    nS = nanmax(sess);
    for iM=1:nM
        for iS=1:nS
            k = mouse==iM & sess==iS;
            k2 = mouseT==iM & sessT==iS;
            casting = find(znC(k) > 0.5);
            iC=1;
            if ~isempty(k)
                while iC<length(casting);
                    if casting(iC)>101;
                        k0 = casting(iC)-100:casting(iC)+100;
                        plot(xT1(k2),yT1(k2),'r.');
                        hold on;
                        plot(x(k),y(k),'.','markersize',1,'color',[0.5 0.5 0.5]);
                        plot(x(k0),y(k0),'k.'); scatter(nx(k0),ny(k0),[],znC(k0),'filled');
                        set(gca,'XLim',[1 1250],'YLim',[1 1000]);
                        colorbar;
                        caxis([-0.01 0.01]);
                        pause;
                        clf;
                        iC=iC+101;
                    else
                        iC=iC+1;
                    end
                end
            end
        end
    end
else
    casting = find(znidphi > 0.5); %#ok<UNRCH>
    iC=1;
    while iC<length(casting);
        if casting(iC)>101;
            k0 = casting(iC)-100:casting(iC)+100;
            plot(yT1,xT1,'r.');
            hold on;
            plot(x0,y0,'.','markersize',1,'color',[0.5 0.5 0.5]);
            plot(x(k0),y(k0),'k.'); scatter(nx(k0),ny(k0),[],znidphi(k0));
            set(gca,'XLim',[1 1250],'YLim',[1 1000]);
            colorbar;
            caxis([0 1]);
            pause;
            clf;
            iC=iC+101;
        else
            iC=iC+1;
        end
    end
end