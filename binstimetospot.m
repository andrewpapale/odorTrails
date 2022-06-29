mint = nanmin((t2s-tff)./initD0); 
maxt = 2.5*nanmean((t2s-tff)./initD0); 
bins = linspace(mint,maxt,20); 
H0 = hist((t2s(conc0==0)-tff(conc0==0))./initD0(conc0==0),bins);
H1 = hist((t2s(conc0==1)-tff(conc0==1))./initD0(conc0==1),bins); 
H2 = hist((t2s(conc0==2)-tff(conc0==2))./initD0(conc0==2),bins); 
H3 = hist((t2s(conc0==3)-tff(conc0==3))./initD0(conc0==3),bins); 
H4 = hist((t2s(conc0==4)-tff(conc0==4))./initD0(conc0==4),bins);


H = histcn([group0,conc0],1:6,0:4,'AccumData',t2s-tff,'fun',@nanmean);
dH = histcn([group0,conc0],1:6,0:4,'AccumData',t2s-tff,'fun',@nanstderr);