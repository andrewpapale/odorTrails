for im=10:20
    for id=1:10
        figure(1); clf;
        subplot(1,3,1);
        plot3(w,y,abs(dphi(:,im,id)),'k'); hold on;
        scatter3(w,y,abs(dphi(:,im,id)),1,abs(dphi(:,im,id)));
        try
            caxis([0 0.5*nanmean(abs(dphi(:,im,id)))]);
        catch
        end
        title('dphi');
        subplot(1,3,2);
        plot3(w,y,abs(dphi1(:,im,id)),'k'); hold on;
        scatter3(w,y,abs(dphi1(:,im,id)),1,abs(dphi1(:,im,id)));
        set(gca,'ZLim',[0 1000]);
        try
            caxis([0 0.5*nanmean(abs(dphi1(:,im,id)))]);
        catch
        end
        xlabel('angular freq');
        ylabel('amplitude');
        zlabel('Curvature');
        title('dphi analytical');
        subplot(1,3,3);
        plot3(w,y,C(:,im,id),'k'); hold on;
        scatter3(w,y,C(:,im,id),1,C(:,im,id));
        set(gca,'ZLim',[0 1000]);
        try
            caxis([0 0.5*nanmean(C(:,im,id))]);
        catch
        end
        title('tortuosity');
        drawnow;
        disp(d(id));
        disp(m(im));
        pause;
    end
end

for id=1:20
    figure(1); clf;
    subplot(1,3,1);
    plot3(w,y,abs(dphi2(:,id)),'k'); hold on;
    scatter3(w,y,abs(dphi2(:,id)),1,abs(dphi2(:,id)));
    try
        caxis([0 0.5*nanmean(abs(dphi2(:,id)))]);
    catch
    end
    title('dphi');
    subplot(1,3,2);
    plot3(w,y,abs(dphi12(:,id)),'k'); hold on;
    scatter3(w,y,abs(dphi12(:,id)),1,abs(dphi12(:,id)));
    try
        caxis([0 0.5*nanmean(abs(dphi12(:,id)))]);
    catch
    end
    xlabel('angular freq');
    ylabel('amplitude');
    zlabel('Curvature');
    title('dphi analytical');
    subplot(1,3,3);
    plot3(w,y,C2(:,id),'k'); hold on;
    scatter3(w,y,C2(:,id),1,C2(:,id));
    try
        caxis([0 0.5*nanmean(C2(:,id))]);
    catch
    end
    title('tortuosity');
    drawnow;
    disp(windowL(id));
    pause;
end


k = abs(y)>6;
figure(1); clf;
for id=1:2:10
    for im=11%1:5:20
        
        subplot(1,3,1);
        plot(w(k),(abs(dphi(k,im,id))),'.'); hold on;
        colorbar;
        set(gca,'YLim',[0 200]);
        set(gca,'XLim',[0 max(w(k))]);
        title('dphi');
        subplot(1,3,2);
        plot(w(k),(abs(dphi1(k,im,id))),'.'); hold on;
        colorbar;
        set(gca,'YLim',[0 200]);
        set(gca,'XLim',[0 max(w(k))]);
        title('dphi analytical');
        subplot(1,3,3);
        plot(w(k),(C(k,im,id)),'.'); hold on;
        colorbar;
        set(gca,'YLim',[0 200]);
        set(gca,'XLim',[0 max(w(k))]);
        title('Tortuosity');
        drawnow;
        disp(d(id));
        xlabel('ang freq','fontsize',18);
        %pause;
    end
end
subplot(1,3,1); lsline;
subplot(1,3,2); lsline;
subplot(1,3,3); lsline;
legend(mat2str(d(1)),mat2str(d(3)),...
    mat2str(d(5)),mat2str(d(7)),mat2str(d(9)));

k = abs(y)>6;
figure(1); clf;
for im=1:5:20
    subplot(1,3,1);
    plot(w(k),(abs(dphi2(k,im))),'.'); hold on;
    colorbar;
    line([0 8],[0 200],'color','k');
    set(gca,'YLim',[0 200]);
    set(gca,'XLim',[0 max(w(k))]);
    title('dphi');
    subplot(1,3,2);
    plot(w(k),(abs(dphi12(k,im))),'.'); hold on;
    colorbar;
    line([0 8],[0 200],'color','k');
    set(gca,'YLim',[0 200]);
    set(gca,'XLim',[0 max(w(k))]);
    title('dphi analytical');
    subplot(1,3,3);
    plot(w(k),(C2(k,im)),'.'); hold on;
    colorbar;
    line([0 8],[0 200],'color','k');
    set(gca,'YLim',[0 200]);
    set(gca,'XLim',[0 max(w(k))]);
    title('Tortuosity');
    drawnow;
    disp(d(id));
end
subplot(1,3,1); lsline;
subplot(1,3,2); lsline;
subplot(1,3,3); lsline;



for id=9%1:10
    for im=11%1:20
        figure(1); clf;
        subplot(1,3,1);
        plot(abs(y),squeeze(abs(dphi(:,im,1:2:10))),'.');
        colorbar;
        lsline;
        line([0 8],[0 200],'color','k');
        set(gca,'YLim',[0 200]);
        %set(gca,'XLim',[0 max(w(k))]);
        title('dphi');
        subplot(1,3,2);
        plot(abs(y),squeeze(abs(dphi1(:,im,1:2:10))),'.');
        colorbar;
        lsline;
        line([0 8],[0 200],'color','k');
        set(gca,'YLim',[0 200]);
        %set(gca,'XLim',[0 max(w(k))]);
        title('dphi analytical');
        subplot(1,3,3);
        plot(abs(y),squeeze(C(:,im,1:2:10)),'.');
        colorbar;
        lsline;
        line([0 8],[0 200],'color','k');
        set(gca,'YLim',[0 200]);
        %set(gca,'XLim',[0 max(w(k))]);
        title('Tortuosity');
        drawnow;
        disp(d(id));
        %pause;
    end
end
legend(mat2str(d(1)),mat2str(d(3)),...
    mat2str(d(5)),...
    mat2str(d(7)),mat2str(d(9)));


subplot(1,3,1);
plot3(w,y,abs(dphi(:,:)));
title('dphi'); subplot(1,3,2);
plot3(w,y,abs(dphi1(:,:)));
title('dphi analytical');
subplot(1,3,3);
plot3(w,y,C(:,:));
title('tortuosity');
pause;

[~,pks] = findpeaks(abs(y));

for id=1:size(dphi,3)
    figure(1); clf;
    plot(w(pks),abs(C(pks,:,id)));
    legend(mat2str(m(1)),mat2str(m(2)),mat2str(m(3)),mat2str(m(4)),mat2str(m(5)),mat2str(m(6)),mat2str(m(7)),mat2str(m(8)),mat2str(m(9)),mat2str(m(10)));
    title(mat2str(d(id)));
    drawnow;
    pause;
end


for id=1:10
    figure(1); clf;
    var = squeeze(abs(dphi1(:,11,id)));
    H = histcn([y',w'],linspace(min(y),max(y),100),linspace(min(w),max(w),100),...
        'AccumData',var,'fun',@nanmean);
    imagesc(H);
    caxis([0 0.5*nanmean(H(:))]);
    pause;
    title(mat2str(d(id)));
%     F = figure(1);
%     saveStr = strcat('foaw-idphi analytical-',mat2str(id),'.jpg');
%     saveas(F,saveStr);
end

iC=1;
ps = postSmoothing;
[~,~,IPana]=InstFreq(1/50,y_dot);
for is=1:1:10
    [~,~,IP0]=InstFreq(1/50,dy0(:,is));
    IP(:,iC) = IP0;
    iC=iC+1;
end
H=[]; for ih=1:1:10; H0 = nansum(IP(:,ih)-IPana'); H = cat(1,H,H0); end;
dH=[]; for ih=1:1:10; dH0 = nanstd(IP(:,ih)-IPana'); dH = cat(1,dH,dH0); end;
bar(ps,H); hold on;
errorbar(ps,H,dH,dH,'.');
