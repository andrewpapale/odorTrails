function compareHists(znC,bait,conc,k,nbins)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

med = nanmedian(znC);
bins = linspace(nanmax(med-20*med,nanmin(znC)),med+20*med,nbins);
bins = unique(bins);

C00 = hist(znC(bait==0 & conc==0 & k),bins);
C01 = hist(znC(bait==0 & conc==1 & k),bins);
C02 = hist(znC(bait==0 & conc==2 & k),bins);
C10 = hist(znC(bait==1 & conc==0 & k),bins);
C11 = hist(znC(bait==1 & conc==1 & k),bins);
C12 = hist(znC(bait==1 & conc==2 & k),bins);

F = figure(1); clf;
subplot(1,2,1);
title('unbaited');
plot(bins,cumsum(C00./nansum(C00))); hold on;
plot(bins,cumsum(C01./nansum(C01)));
plot(bins,cumsum(C02./nansum(C02)));
legend('0.1%','1%','2%');

subplot(1,2,2);
title('baited');
plot(bins,cumsum(C10./nansum(C10))); hold on;
plot(bins,cumsum(C11./nansum(C11)));
plot(bins,cumsum(C12./nansum(C12)));
legend('0.1%','1%','2%');

F = figure(2); clf;
subplot(1,2,1);
title('unbaited');
H = cat(1,C00./nansum(C00),C01./nansum(C01),C02./nansum(C02));
bar(bins,H'); hold on;
legend('0.1%','1%','2%');

subplot(1,2,2);
title('baited');
H = cat(1,C10./nansum(C10),C11./nansum(C11),C12./nansum(C12));
bar(bins,H'); hold on;
legend('0.1%','1%','2%');

F = figure(3);
H = histcn([conc(k),bait(k)],0:2,0:1,'AccumData',znC(k),'fun',@nanmean);
bar(H);

F = figure(4);
boxplot(znC(k),{conc(k),bait(k)},'notch','on','color','bybyby');
set(gca,'YLim',[nanmedian(znC(k))-10 nanmedian(znC(k))+10]);
end

