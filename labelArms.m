function labelArms

ftrail = dir('*-Odor-Trail.mat');
fcp = dir('*-ConnectingPoint.mat');

assert(length(ftrail)==length(fcp),'file mismatch');

nF = length(fcp);
for iF=1:nF
    
    load(ftrail(iF).name);
    load(fcp(iF).name);
    
    xT0 = yT;
    yT0 = xT;
    
    outside1 = zeros(size(data));
    xTO = xT(outside==1);
    yTO = yT(outside==1);
    for iT=1:length(xTO)
        outside1(xTO(iT),yTO(iT))=1;
    end
    
    CC = bwconncomp(outside1);
    
    disp(CC.NumObjects);
    for iR=1:CC.NumObjects
        if length(CC.PixelIdxList{iR})>100
            [i,j] = ind2sub(size(data),CC.PixelIdxList{iR});
            keepV = zeros(length(xT0),1);
            for iP=1:length(xT0)
                for iQ=1:length(j)
                    if xT0(iP)==j(iQ) & yT0(iP)==i(iQ); %#ok<AND2>
                        keepV(iP)=1;
                    end
                end
            end
            
            figure(1); clf;
            plot(xT0,yT0,'r.'); hold on;
            plot(xT0(keepV==1),yT0(keepV==1),'g.');
            quiver(xC0{1},yC0{1},xC-xC0{1},yC-yC0{1},'b');
            legend('Trail','Arm to Label','Arm 1');
            
            iD = input('Enter 1,2,or 3 :');
            
            arm{iD} = keepV;
        end
    end
    
    saveStr = ftrail(iF).name;
    disp(saveStr);
    save(saveStr,'xC','yC','outside','inside','vec','xC0','yC0','arm');
    
end
