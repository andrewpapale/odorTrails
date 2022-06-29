function tempVideoReader
[handles.FileName,pathName,~] = uigetfile;

extStr = strsplit(handles.FileName,'.');
assert(strcmp(extStr{2},'avi'),'file must be an avi video file');
v = VideoReader(handles.FileName);
video = read(v,[1 250]); %#ok<VIDREAD>

% get matching position file
cd(strcat(pathName,'\positions')); % this is also where file will be saved
tempStr = strsplit(extStr{1},'_');



flag = false;
while ~flag
    iFrame = input('Enter frame: ');
    
    figure(1); clf;
    imagesc(squeeze(video(:,:,:,iFrame)));
    hold on;
    axis xy
    x = position_results.mouseCOM(iFrame,1);
    y = position_results.mouseCOM(iFrame,2);
    nx = position_results.nosePOS(iFrame,1);
    ny = position_results.nosePOS(iFrame,2);
    plot(x,y,'b.','markersize',20);
    plot(nx,ny,'r.','markersize',20);
    title(mat2str(iFrame));
    
    okStr = input('OK? ','s');
    if strcmp(okStr,'y')
        flag = true;
        tempStr = strsplit(posFile.name,'-');
        timeStr = tempStr{4};
        disp(timeStr);
        saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-','StartFrame.mat');
        startFrame = iFrame;
        save(saveStr,'startFrame');
    end
    
end


end