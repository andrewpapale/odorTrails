% determine which sessions are excluded
cd('C:\Users\papalea\Documents\Data');
[num,text] = xlsread('KP_AL_CombinedLog.xlsx',1,'A1:H245','basic');
count = 0;
for iC=1:5;
    for iT=1:231;
        currcell = text(iT,iC);
        if ~strcmp(currcell,'');
            count = count+1;
        end
    end
end

cd('C:\Users\papalea\Documents\Data\Spot\positions');
fp = dir('*positions.mat');


% mouseStr
iC = 1;
for iT=1:245
    mouse = mat2str(num(iC,2));
    date0 = mat2str(num(iC,1));
    mouseStr =  strcat(mouse,'_','20',date0(1:2),'-',date0(3:4),'-',date0(5:6));
    
    % update mouseStr?
    for iD=1:5
        if ~strcmp(text(iT,iD),'')
            % don't update mouseStr
        else
            % update mouseStr
            iC = iC+1;
            mouse = mat2str(num(iC,2));
            date0 = mat2str(num(iC,1));
            mouseStr =  strcat(mouse,'_','20',date0(1:2),'-',date0(3:4),'-',date0(5:6));
        end
    end
    
    fp0 = fp(iT).name;
    fp0 = fp0(1:17);
    disp(iT);
    if strcmp(fp0,mouseStr);
    else
        disp(mouseStr);
        disp(fp(iT).name);
        fprintf('trial = %d \n',iD);
        pause;
        cd('C:\Users\papalea\Documents\Data');
        [num,text] = xlsread('KP_AL_CombinedLog.xlsx',1,'A1:H245','basic');
    end
    
end