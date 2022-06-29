figure; clf;

k = C==0; 
pR = R(:,k);
N = nansum(pR,1);
mpR = nanmean(pR./(repmat(N,[size(pR,1),1]).*area),2);
dpR = nanstderr(pR./(repmat(N,[size(pR,1),1]).*area),[],2);
lineProps.col = {'c'};
mseb(radii,log10(mpR'),dpR'./(mpR'*log(10)),lineProps);

k = C==1; 
pR = R(:,k);
N = nansum(pR,1);
mpR = nanmean(pR./(repmat(N,[size(pR,1),1]).*area),2);
dpR = nanstderr(pR./(repmat(N,[size(pR,1),1]).*area),[],2);
lineProps.col = {'b'};
mseb(radii,log10(mpR'),dpR'./(mpR'*log(10)),lineProps);

k = C==2; 
pR = R(:,k);
N = nansum(pR,1);
mpR = nanmean(pR./(repmat(N,[size(pR,1),1]).*area),2);
dpR = nanstderr(pR./(repmat(N,[size(pR,1),1]).*area),[],2);
lineProps.col = {'g'};
mseb(radii,log10(mpR'),dpR'./(mpR'*log(10)),lineProps);

k = C==3; 
pR = R(:,k);
N = nansum(pR,1);
mpR = nanmean(pR./(repmat(N,[size(pR,1),1]).*area),2);
dpR = nanstderr(pR./(repmat(N,[size(pR,1),1]).*area),[],2);
lineProps.col = {'m'};
mseb(radii,log10(mpR'),dpR'./(mpR'*log(10)),lineProps);

k = C==4; 
pR = R(:,k);
N = nansum(pR,1);
mpR = nanmean(pR./(repmat(N,[size(pR,1),1]).*area),2);
dpR = nanstderr(pR./(repmat(N,[size(pR,1),1]).*area),[],2);
lineProps.col = {'r'};
mseb(radii,log10(mpR'),dpR'./(mpR'*log(10)),lineProps);

