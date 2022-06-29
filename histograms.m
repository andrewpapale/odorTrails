k = ~edge & dnT<=15 & towards;
bins = linspace(2.5,100,10);

H = [];
H1 = [];
for iC=bins;
    k0 = k & nC < 10;
    H0 = histcn([bait(k0),conc(k0)],0:1,0:2,'AccumData',znidphi(k0),'fun',@nanmean);
    H = cat(3,H,H0);
end
for iC=bins; 
    k0 = k & nC < 10;
    H0 = histcn([bait(k0),conc(k0)],0:1,0:2,'AccumData',znC(k0),'fun',@nanmean);
    H1 = cat(3,H1,H0);
end
F1 = figure(1);
subplot(1,2,1); 
imagesc(1:3,bins,squeeze(H(1,:,:))'); 
colorbar; set(gca,'YTick',bins); 
subplot(1,2,2); 
imagesc(1:3,bins,squeeze(H1(1,:,:))'); 
set(gca,'YTick',bins); 
colorbar;

k = ~edge;

H = [];
H1 = [];
for iC=bins;
    k0 = k & nC<10;
    H0 = histcn([conc(k0),dnT(k0)],0:2,cat(2,quantile(dnT(k0),49),Inf),'AccumData',znidphi(k0),'fun',@nanmean);
    H = cat(3,H,H0);
end

for iC=bins; 
    k0 = k & nC < 10;
    H0 = histcn([conc(k0),dnT(k0)],0:2,cat(2,quantile(dnT(k0),49),Inf),'AccumData',znC(k0),'fun',@nanmean);
    H1 = cat(3,H1,H0);
end

F2 = figure(2);
subplot(3,2,1);
imagesc(squeeze(H(1,:,:))');
colorbar; set(gca,'YTick',bins); 
subplot(3,2,2);
imagesc(squeeze(H(2,:,:))');
colorbar; set(gca,'YTick',bins); 
subplot(3,2,3);
imagesc(squeeze(H(3,:,:))');
colorbar; set(gca,'YTick',bins); 

F2 = figure(2);
subplot(3,2,4);
imagesc(squeeze(H1(1,:,:))');
colorbar; set(gca,'YTick',bins); 
subplot(3,2,5);
imagesc(squeeze(H1(2,:,:))');
colorbar; set(gca,'YTick',bins); 
subplot(3,2,6);
imagesc(squeeze(H1(3,:,:))');
colorbar; set(gca,'YTick',bins); 

k = ~edge & towards;
k0 = k;
H = histcn([conc(k0),dnT(k0)],0:2,cat(2,quantile(dnT(k0),49),Inf),'AccumData',znidphi(k0),'fun',@nanmean);
k0 = k & nC<10;
H1 = histcn([conc(k0),dnT(k0)],0:2,cat(2,quantile(dnT(k0),49),Inf),'AccumData',znC(k0),'fun',@nanmean);

F3 = figure(3);
subplot(1,2,1);
plot(quantile(dnT(k0),49),H');
set(gca,'XLim',[0 30]);
subplot(1,2,2);
plot(quantile(dnT(k0),49),H1');
set(gca,'XLim',[0 30]);

k = ~edge & towards & dnT<=15;
F4 = figure(4);
k0 = k;
H = histcn([bait(k0),conc(k0)],0:1,0:2,'AccumData',znidphi(k0),'fun',@nanmean);
k0 = k;
H1 = histcn([bait(k0),conc(k0)],0:1,0:2,'AccumData',znC(k0),'fun',@nanmean);
subplot(1,2,1); bar(H');
subplot(1,2,2); bar(H1');

znC = nan(size(nC));
k0 = nC < 0.0001;
for iM=1:max(mouse);
    for iS=1:max(sess);
        k1 = mouse==iM&sess==iS;
        if sum(k1)>0
            znC(k1 & k0)=nanzscore(nC(k1 & k0));
        end
    end
end

