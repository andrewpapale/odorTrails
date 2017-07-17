function getTrailComponents(nSegments)
% 2017-06-22 AndyP
% getTrailComponents
% function to mark line segments and connected portion of trail

[handles.FileName,~,~] = uigetfile;
load(handles.FileName);

figure(1); clf;
imagesc(data); axis xy;
hold on;

radius = 100;
boxSize = 100;
nSegments = 3;

okflag = 0;
while ~okflag
    fprintf('Click on connecting point \n');
    [xC,yC] = ginput(1);
    if ~isempty(xC)
        % find closest endpoint to click
        P1 = plot(xC,yC,'gx','markersize',10); hold on;
        P2 = plot(xC,yC,'go','markersize',10); hold on;
    else
        disp('Try again');
    end
    % draw analysis circle around connecting point
    
    H = circle([xC,yC],radius,1000,'r-');
    
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
            delete(P2);
            delete(H);
        else
            flag3 = 0;
        end
    end
    
    
end

[yD,xD] = find(data==1);
outside = sqrt((xD-xC).^2+(yD-yC).^2)>radius;
inside = sqrt((xD-xC).^2+(yD-yC).^2)<=radius;

P1 = plot(yT(outside),xT(outside),'g.');
P2 = plot(yT(inside),xT(inside),'r.');

pause;
delete(P1);
delete(P2);


for iS=1:nSegments
    figure(1); clf;
    imagesc(data); axis xy;
    hold on;
    okflag = 0;
    while ~okflag
        fprintf('Click on line segments \n');
        [xC2,yC2] = ginput(2);
        if ~isempty(xC)
            % find closest endpoint to click
            P1 = plot(xC2(1),yC2(1),'gx','markersize',10); hold on;
            P2 = plot(xC2(1),yC2(1),'go','markersize',10); hold on;
            P3 = plot(xC2(2),yC2(2),'gx','markersize',10); hold on;
            P4 = plot(xC2(2),yC2(2),'go','markersize',10); hold on;
        else
            disp('Try again');
        end
        
        m = (yC2(2) - yC2(1)) / (xC2(2) - xC2(1));
        if m<1 && m>0
            m=0.1;
        elseif m>-1 && m<0
            m=-0.1;
        end
        mp = -1/m;
        result = computeline([xC2(1),yC2(1)],mp,[xC2(1)-boxSize xC2(1)+boxSize]);
        for iR=1:length(result)
            result1(iR,1:2) = result{iR};
        end
        result = computeline([xC2(2),yC2(2)],mp,[xC2(2)-boxSize xC2(2)+boxSize]);
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
        k = tempX<=1 | tempY<=1 | tempX>size(data,2) | tempY>size(data,1) | isnan(tempX) | isnan(tempY);
        tempX(k)=[];
        tempY(k)=[];
        for iL=1:length(tempX)
            tempsegment(round(tempY(iL)),round(tempX(iL)))=1;
        end
        se = strel('line',11,90);
        tempsegment = imdilate(tempsegment,se);
        tempsegment = bwmorph(tempsegment,'majority');
        [yS0,xS0]=find(tempsegment==1);
        [ySI0,xSI0] = find(data==1 & tempsegment==1);
        
        keyboard;
        
        P5 = plot(xS0,yS0,'gx');
        P6 = plot(xSI0,ySI0,'r.');
        
        % compute vector (pointing towards connecting endpoint)
        x1 = xC-xC2;
        y1 = yC-yC2;
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
                delete(P5);
                delete(P6);
            else
                flag3 = 0;
            end
        end
        
    end
end

tempStr = strsplit(handles.FileName,'-');
saveStr = strcat(tempStr{1},'-',tempStr{2},'-',tempStr{3},'-',tempStr{4},'-',tempStr{5},'-','ConnectingPoint.mat');
save(saveStr,'xC','yC','outside','inside','xS','yS','segment','xSI','ySI','vec','xC0','yC0');


end


