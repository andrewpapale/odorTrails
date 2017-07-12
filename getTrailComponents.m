function getTrailComponents(nSegments)
% 2017-06-22 AndyP
% getTrailComponents
% function to mark line segments and connected portion of trail

[handles.FileName,~,~] = uigetfile;
load(handles.FileName);

figure(1); clf;
imagesc(data); axis xy;
hold on;

boxSize = 100;

okflag = 0;
while ~okflag
    fprintf('Click on connecting point \n');
    [xC,yC] = ginput(1);
    if ~isempty(xC)
        % find closest endpoint to click
        P1 = plot(xC,yC,'gx','markersize',10); hold on;
        P2 = plot(xC,yC,'go','markersize',10); hold on;
        X= round(xC);
        Y= round(yC);
    else
        disp('Try again');
    end
    % draw analysis rectangles around each segment and box around connected
    % endpoint
    box=zeros(size(data));
    box(Y-boxSize:Y+boxSize,X-boxSize:X+boxSize)=1;
    
    [yB,xB]=find(box==1);
    P1 = plot(xB,yB,'rx');
    
    
    % ask user if this is OK
    flag3 = 0;
    while ~flag3
        okStr = input('OK? (Enter "y" or "n")  ','s');
        if strcmp(okStr, 'y');
            okflag = 1;
            flag3 = 1;
        elseif strcmp(okStr, 'n');
            okflag = 0;
            flag3 = 1;
            delete(P1);
        else
            flag3 = 0;
        end
    end
    
    
end

[yTO,xTO] = find(data==1 & box~=1);
[yTI,xTI] = find(data==1 & box==1);


for iS=1:nSegments
    figure(1); clf;
    imagesc(data); axis xy;
    hold on;
    okflag = 0;
    while ~okflag
        fprintf('Click on line segments \n');
        [xC,yC] = ginput(2);
        if ~isempty(xC)
            % find closest endpoint to click
            P1 = plot(xC,yC,'gx','markersize',10); hold on;
            P2 = plot(xC,yC,'go','markersize',10); hold on;
        else
            disp('Try again');
        end
        
        m = (yC(2) - yC(1)) / (xC(2) - xC(1));
        mp = -1/m;
        result = computeline([xC(1),yC(1)],mp,[xC(1)-boxSize xC(1)+boxSize]);
        for iR=1:length(result)
            result1(iR,1:2) = result{iR};
        end
        result = computeline([xC(2),yC(2)],mp,[xC(2)-boxSize xC(2)+boxSize]);
        for iR=1:length(result)
            result2(iR,1:2) = result{iR};
        end
        for iR=1:length(result)
            line{iR} = linepts([result1(iR,1),result1(iR,2)],[result2(iR,1),result2(iR,2)]);
        end
        tempX = [];
        tempY = [];
        for iR=1:length(result)
            templine = line{iR};
            nL=length(templine);
            for iL=1:nL
                tempX = cat(1,tempX,templine{iL}(1,1));
                tempY = cat(1,tempY,templine{iL}(1,2));
            end
        end
        tempsegment = zeros(size(data));
        k = tempX<=0 | tempY<=0;
        tempX(k)=[];
        tempY(k)=[];
        for iL=1:length(tempX)
            tempsegment(round(tempY(iL)),round(tempX(iL)))=1;
        end
        [yS0,xS0]=find(tempsegment==1);
        [ySI0,xSI0] = find(data==1 & tempsegment==1);
        
        
        P3 = plot(xS0,yS0,'gx');
        P4 = plot(xSI0,ySI0,'r.');
        
        % compute vector (pointing towards connecting endpoint)
        x1 = xC(1)-xC(2);
        y1 = yC(1)-yC(2);
        mag = sqrt(x1.^2+y1.^2);
        
        
        % ask user if this is OK
        flag3 = 0;
        while ~flag3
            okStr = input('OK? (Enter "y" or "n")  ','s');
            if strcmp(okStr, 'y');
                okflag = 1;
                flag3 = 1;
                xS{iS} = xS0;
                yS{iS} = yS0;
                segment{iS} = tempsegment;
                xSI{iS} = xSI0;
                ySI{iS} = ySI0;
                vec.x0{iS} = x1;
                vec.y0{iS} = y1;
                vec.mag{iS} = mag;
                xC0{iS} = xC;
                yC0{iS} = yC;
            elseif strcmp(okStr, 'n');
                okflag = 0;
                flag3 = 1;
                delete(P1);
                delete(P2);
                delete(P3);
                delete(P4);
            else
                flag3 = 0;
            end
        end
        
    end
end

tempStr = strsplit(handles.FileName,'-');
saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-',tempStr{5},'-','ConnectingPoint.mat');
save(saveStr,'box','xB','yB','xTI','yTI','xTO','yTO','X','Y','xS','yS','segment','xSI','ySI','vec','xC0','yC0');


end


