% 2019-11-11 AndyP
% Tara_CutOut2

% open image
[FileName,~,~] = uigetfile;
I = imread(FileName);
F = figure(1); clf;
imagesc(I);
hold on;

% user inputs layer that they are counting
name = input('Layer 1-4, Layer 5-6 ','s');

inStr = input('Do you want to load a previously drawn line?   ','s');
if strcmp(inStr,'yes')
    fprintf('load previously drawn line \n');
    [fn1,~,~] = uigetfile;
    load(fn1);
    plot(roi.Position(:,1),roi.Position(:,2),'r.-','markersize',30);
else
end



roi = images.roi.Polyline(gca);
draw(roi);





% fprintf('open .mat file for PV \n');
% uiopen;

% [inPV,onPV] = inpolygon(cell_pos1(:,1),cell_pos1(:,2),roi.Position(:,1),roi.Position(:,2));
% nPV = sum(inPV);
% 
% fprintf('open .mat file for VVA \n');
% uiopen;
% 
% [inVVA,onVVA] = inpolygon(cell_pos1(:,1),cell_pos1(:,2),roi.Position(:,1),roi.Position(:,2));
% nVVA = sum(inVVA);
% 
% dateStr = strsplit(datestr(now),' ');
% tempStr = strsplit(FileName,'_');
% saveStr = strcat(tempStr{1},'_',name,'_','cellcounts','_',dateStr{1},'.mat');
% save(saveStr,'roi','inPV','inVVA','onPV','onVVA','nPV','nVVA');

dateStr = strsplit(datestr(now),' ');
tempStr = strsplit(FileName,'_');
tempStr = tempStr{1};
tempStr = strsplit(tempStr,'-');
tempStr = strcat(tempStr{1},'-',tempStr{2});
saveStr = strcat(tempStr,'_',name,'_','CutOut','_',dateStr{1},'.mat');
save(saveStr,'roi');



