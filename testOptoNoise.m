figure(1); clf;
ix = interp1(time0,frame0,timeOut0(AmpOutNoise==5),'nearest')+1;
ix1 = interp1(time0,frame0,timeOut0(AmpOutNoise==0),'nearest')+1;
f2s = interp1(time0,frame0,foundSpotT,'nearest');
f2f = interp1(time0,frame0,trialStartT,'nearest');

subplot(8,8,[5:8,13:16,21:24]);
scatter(x0,y0,10,time0,'filled');
hold on;
plot(x0(ix),y0(ix),'rx','linewidth',2,'markersize',10);
plot(x0(ix1),y0(ix1),'bx','linewidth',2,'markersize',10);
plot(x0(f2s),y0(f2s),'co','linewidth',2,'markersize',10);
plot(x0(f2f),y0(f2f),'mo','linewidth',2,'markersize',10);
viscircles([cx,cy],2*5.174,'color','r');
set(gca,'fontsize',18);
axis off;
XL = get(gca,'XLim');
YL = get(gca,'YLim');


subplot(8,8,[1:4,9:12,17:20]);
plot(time0,y0,'k.');
hold on;
plot(time0(ix),y0(ix),'rx','linewidth',2,'markersize',10);
plot(time0(ix1),y0(ix1),'bx','linewidth',2,'markersize',10);
line([time0(f2s) time0(f2s)],[0 600],'color','c');
line([time0(f2f) time0(f2f)],[0 600],'color','m');
xlabel('time (s)');
set(gca,'fontsize',18);
%h = gca; h.YAxis.Visible = 'off';
set(gca,'YLim',YL);

subplot(8,8,[29:32,37:40,45:48,53:56]);
plot(x0,time0,'k.');
hold on;
plot(x0(ix),time0(ix),'rx','linewidth',2,'markersize',10);
plot(x0(ix1),time0(ix1),'bx','linewidth',2,'markersize',10);
line([0 600],[time0(f2s) time0(f2s)],'color','c');
line([0 600],[time0(f2f) time0(f2f)],'color','m');
ylabel('time (s)');
set(gca,'fontsize',18);
%h = gca; h.XAxis.Visible = 'off';
set(gca,'XLim',XL);