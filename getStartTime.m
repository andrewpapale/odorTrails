function getStartTime

homeDir = cd;

movieFiles = dir('*.avi');
nM = length(movieFiles);

for iM=1:nM
    
    
    
    extStr = strsplit(movieFiles(iM).name,'.');
    
    
    % get matching position file
    cd(strcat(homeDir,'\positions')); % this is also where file will be saved
    
    posFile = dir(strcat(extStr{1},'*_positions.mat'));
    if ~isempty(posFile)
        tempStr = strsplit(posFile.name,'-');
        saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-','StartFrame.mat');
        alreadyDone = dir(saveStr);
    else
        alreadyDone.name = '000';
    end
    
    if ~isempty(alreadyDone)
        if strcmp(alreadyDone.name,'000');
            disp(extStr{1});
        end
    end
    
    if isempty(alreadyDone)
        
        load(posFile.name,'position_results');
        
        cd(homeDir);
        
        % load video file
        v = VideoReader(movieFiles(iM).name);
        video = read(v,[1 500]); %#ok<VIDREAD>
        
        flag = false;
        while ~flag
            flag2 = false;
            while ~flag2
                iFrame = input('Enter frame: ');
                if ~isempty(iFrame)
                    if iFrame >=1 && iFrame <=500
                        flag2=true;
                    end
                end
            end
            
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
                saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-','StartFrame.mat');
                startFrame = iFrame;
                save(saveStr,'startFrame');
                disp(saveStr);
            end
            
        end
    end
    
    clear alreadyDone
end