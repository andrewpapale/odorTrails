pos = dir('*_positions.mat');

nP = length(pos);

for iP=1:nP
    tempStr = strsplit(pos(iP).name,'-');
    timeStr = tempStr{4};
    disp(timeStr);
    saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-','StartFrame.mat');
    int = input('Enter Start Time:   ');
    startFrame = int;
    save(saveStr,'startFrame');
end