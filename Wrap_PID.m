function [meanV,medianV,stdV,TablePosition,PotBias,fracPos,Condition,trace] = Wrap_PID
% 2019-06-01 AndyP

[num,txt] = xlsread('2019-06-11-PID.xlsx','A2:G38');
% 
D = dir('*.mat');
nD = length(D);
% nD = length(D)-5;
% 
TablePosition = num(1:nD,2);
PotBias = num(1:nD,5);
% TablePosition = [];
% PotBias = [];

meanV = [];
medianV = [];
stdV = [];
fracPos = [];
trace = [];
for iD=1:nD
    load(D(iD).name,'data');
    meanV = cat(1,meanV,nanmean(data));
    medianV = cat(1,medianV,nanmedian(data));
    stdV = cat(1,stdV,nanstd(data));
    fracPos = cat(1,fracPos,sum(data-nanmean(data)>0)./length(data));
    
    trace = cat(2,trace,data);
    
    if strcmp(txt{iD,5},'E')
        Condition(iD)=1;
    else strcmp(txt{iD,5},'B')
        Condition(iD)=0;
    end
    
end

