% analyze PID script
% 2019-07-12 AndyP

D = dir('*.mat');
nD = length(D);
nP = [];
dist = [];
dtP = [];
pP = [];
mV = [];
meanVar = [];
I = [];
distMatrix = [15 14	13 12 11 10 9 8 7 6 5 4 3 2.5 2 1.5 1 0.5];
%distMatrix = [0	26	27	28	36	35	34	33	25	17	18	19	20	42	0	43	44	50	51	52	60	67	66	65	11	12	10	0	9	8	2	3	4	5	6	14	15	16	24	23	22	30	31	32	38	39	40	46	47	48	49	54	55	56	0	57	58	62	63	64	65	73	72	71	70];
cf = 0.005;
d1 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',1000);
d2 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',70,'HalfPowerFrequency2',74, ...
               'DesignMethod','butter','SampleRate',1000);
s0 = [];           
f0 = [];         
D0 = [];
for iD=1:nD
    if distMatrix(iD)~=0
        load(D(iD).name,'timeStamps','data');
        s1 = movmean(filtfilt(d2,filtfilt(d1,data-nanmin(data))),[200,20]);
        [data1, f, cost] = beads(s1, 1, 0.001, 3, 0.1*cf, 2*cf, 3*cf);
        data1 = data1(587:29332);
        f0 = cat(1,f0,(587:29332)');
        D0 = cat(1,D0,data1);
        s0 = cat(1,s0,repmat(distMatrix(iD),[28745,1]));
        [pks,locs,w,p] = findpeaks(data1,'MinPeakProminence',0.03);
        I = cat(1,I,nanmean(data1>cf));
        nP = cat(1,nP,pks);
        dist = cat(1,dist,repmat(distMatrix(iD),[length(pks),1]));
        dtP = cat(1,dtP,nan,diff(timeStamps(locs)));
        pP = cat(1,pP,p);
        mV = cat(1,mV,nanmean(data1));
        meanVar = cat(1,meanVar,nanmean(var(data1)));
    end
end


% FigureS#A
xl = sort(distMatrix,'ascend');
F = figure(1); clf;
for iD=1:18 
    subplot(1,18,iD); 
    k = s0==xl(iD); 
    plot(f0(k==1),D0(k==1),'k.'); 
    set(gca,'YLim',[0 1]); 
    %title(sprintf('%0.1fs',xl(iD)),'Units','normalized','Position',[0.5 -0.1 0],'fontsize',30); 
    axis off; 
end

% Figure S#B
H = histcn(dist,sort(distMatrix),'AccumData',nP,'fun',@nanmean);
plot(sort(distMatrix),H,'.','markersize',30,'color','k')
hold on;
plot(distMatrix,I,'.','markersize',30,'color',[0.5 0.8 0.1]);
%legend('Mean Peak Height','Intermittency')
xlabel('distance from source (cm)');
ylabel('%');

% Figure S#C
plot((distMatrix),log10(mV),'k.','markersize',30);
ylabel('log_{10} mean output (AU)');
set(gca,'fontsize',18);
xlabel('distance from source (cm)');

