dir = cd;


% load all data
cd(strcat(dir,'\MS on setup table (10s)'));
D = dir('*.mat');
group = [-1 0 1 2 3 4 5 6 7 8 9 10 -1 10 9 8 7 6 5 4 3 2 1 0]; data0 = []; group0 = []; ndata = []; for iD=1:length(D); load(D(iD).name); data0 = cat(1,data0,data); sdata = sort(data,'ascend'); baseline = nanmean(sdata(1:30)); ndata = cat(1,ndata,data-baseline); group0 = cat(1,group0,repmat(group(iD),[length(data),1])); end;

cd(strcat(dir,'\MS towards shelves'));
D = dir('*.mat');
group = [-1 0 1 2 3 4 5 6 7 8 9 10]; for iD=1:length(D); load(D(iD).name); data0 = cat(1,data0,data); sdata = sort(data,'ascend'); baseline = nanmean(sdata(1:30)); ndata = cat(1,ndata,data-baseline); group0 = cat(1,group0,repmat(group(iD),[length(data),1])); end;

cd(strcat(dir,'MS towards sink'));
D = dir('*.mat');
group = [-1 0 1 2 3 4 5 6 7 8 9 10 -1]; for iD=1:length(D); load(D(iD).name); data0 = cat(1,data0,data); sdata = sort(data,'ascend'); baseline = nanmean(sdata(1:30)); ndata = cat(1,ndata,data-baseline); group0 = cat(1,group0,repmat(group(iD),[length(data),1])); end;

[pks,locs,w,p] = findpeaks(ndata,'MinPeakProminence',0.25); % find peaks for normalized data