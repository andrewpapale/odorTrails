bin = quantile(dnT,10);

nXY = ceil(sqrt(length(bin)));
figure(1); clf;
for iT=1:length(bin)-1
    k = ~edge & bin(iT) <= dnT & dnT < bin(iT+1);
    subplot(nXY,nXY,iT);
    circ_plot(spotang(k),'hist',[],100,true,true);
    title(mat2str(round(bin(iT))));
end
  