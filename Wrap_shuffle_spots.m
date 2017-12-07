% Wrap_shuffle_spots

tic;

nshuf = 100;
nsess = 915;

nM = max(mouse);

xT = [];
yT = [];

disp('gathering spots...');
for iM=1:nM
    km = mouse==iM;
    nS = nanmax(sess(km));
    for iS=1:nS
        
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            kT = mouseT==iM & sessT==iS;
            
            % get all spots
            xT = cat(1,xT,nanmedian(xT1(kT)));
            yT = cat(1,yT,nanmedian(yT1(kT)));
            
        end
    end
   
end

initD = nan(915,1000);
t2s = nan(915,1000);
spotfound = nan(915,1000);
dwellT = nan(915,1000);
quadrant = nan(915,1000);
t2s2 = nan(915,1000);
t2s3 = nan(915,1000);
initO = nan(915,1000);
for ishuf = 1:nshuf
    
    % shuffle spots
    rng('shuffle');
    idx = randperm(length(xT));
    xT2 = xT(idx);
    yT2 = yT(idx);
    iC = 1;
    
    disp('computing distances...');
    for iM=1:nM
        km = mouse==iM;
        nS = nanmax(sess(km));
        for iS=1:nS
            
            k = mouse==iM & sess==iS;
            
            if sum(k)>0
                
                
                nx0 = nx(k);
                ny0 = ny(k);
                
                xT0 = xT2(iC);
                yT0 = yT2(iC);
                iC = iC+1;
                
                % calculate distance from trail
                dnT0 = sqrt((nx0-xT0).^2+(ny0-yT0).^2);
                
                shufdnT = cat(1,shufdnT,dnT0);
            end
        end
        
    end
    
    [initD0,t2s0,spotfound0,dwellT0,quadrant0,t2s20,t2s30,initO0]=...
        getTime2Spot2(mouse,sess,shufdnT,x,y,nx,ny,frame,xT2,yT2);
    
    k = ~isnan(t2s0);
    
    initD(:,ishuf) = initD0(k);
    t2s(:,ishuf) = t2s0(k);
    spotfound(:,ishuf) = spotfound0(k);
    dwellT(:,ishuf) = dwellT0(k);
    quadrant(:,ishuf) = quadrant0(k);
    t2s2(:,ishuf) = t2s2(k);
    t2s3(:,ishuf) = t2s3(k);
    initO(:,ishuf) = initO0(k);
    
    fprintf('shuffle %d/%d complete \n',ishuf,nshuf);
    toc;
    
end
