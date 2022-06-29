% zdphi analysis figure

k = ~edge & ~kstop;
bins = quantile(dnT(k),30);
bins1 = quantile(dnT(k),29);
H1 = histcn(dnT(k),bins,'AccumData',znidphi1(k),'fun',@nanmean);
H1(H1==0)=nan;
k = ~edge & ~kstop & t2s1 & spotfound1;
H2 = histcn(dnT(k),bins,'AccumData',znidphi1(k),'fun',@nanmean);
H2(H2==0)=nan;
k = ~edge & ~kstop & ~t2s1 & spotfound1;
H3 = histcn(dnT(k),bins,'AccumData',znidphi1(k),'fun',@nanmean);
H3(H3==0)=nan;
k = ~edge & ~kstop & ~t2s1 & ~spotfound1;
H4 = histcn(dnT(k),bins,'AccumData',znidphi1(k),'fun',@nanmean);
H4(H4==0)=nan;

figure(1); clf; hold on;
plot(bins1,H1);
plot(bins1,H2);
plot(bins1,H3);
plot(bins1,H4);
legend('All','before finding spot','after finding spot','did not find spot');