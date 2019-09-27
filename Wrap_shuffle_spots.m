% Wrap_shuffle_spots

tic;

nshuf = 100;


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

rT0 = xT<44.8 | xT>1280-44.8 | yT<44.8 | yT>1024-44.8;
xT(rT0)=[];
yT(rT0)=[];

nsess = length(xT);

sinitD = nan(nsess,nshuf);
st2s = nan(nsess,nshuf);
sspotfound = nan(nsess,nshuf);
sdwellT = nan(nsess,nshuf);
squadrant = nan(nsess,nshuf);
sinitO = nan(nsess,nshuf);
sinitOb = nan(nsess,nshuf);
sspotring = nan(nsess,nshuf);
ssess = nan(nsess,nshuf);
smouse = nan(nsess,nshuf);
sdtraveled = nan(nsess,nshuf);

for ishuf = 1:nshuf
    
    % shuffle spots
    rng('shuffle');
    idx = randperm(length(xT));

    % correct idx for if spot is the same
    realidx = 1:nsess;
        moveidx = [];
        for iT=1:nsess
            if idx(iT)==realidx(iT)
                moveidx = cat(1,moveidx,idx(iT));
            end
        end
        if length(moveidx)>1
            for iT=1:length(moveidx)-1
                idx(moveidx(iT))=moveidx(iT+1);
            end
            idx(moveidx(end))=moveidx(1);
        elseif length(moveidx)==1
            flag = true;
            while flag
                newidx = randi(nsess);
                if newidx~=moveidx && idx(newidx)~=idx(moveidx)
                    flag = false;
                    tempidx = idx(moveidx);
                    idx(moveidx)=idx(newidx);
                    idx(newidx)=tempidx;
                end
            end
        end
        
        assert(~any(realidx==idx),'error spots not shuffled!');

    
    xT2 = xT(idx);
    yT2 = yT(idx);
    iC = 1;
    
    disp('computing distances...');
    shufdnT = [];
    nx2 = [];
    ny2 = [];
    frame2 = [];
    mouse2 = [];
    sess2 = [];
    x2 = [];
    y2 = [];
    for iM=1:nM
        for iS=1:nS
            
            k = mouse==iM & sess==iS;
            kT = mouseT==iM & sessT==iS;
            
            tempx = nanmedian(xT1(kT));
            tempy = nanmedian(yT1(kT));
            
            if sum(k)>0 && ~(tempx<44.8 | tempx>1280-44.8 | tempy<44.8 | tempy>1024-44.8) %#ok<OR2>
                
                mouse2 = cat(1,mouse2,mouse(k));
                sess2 = cat(1,sess2,sess(k));
                
                x2 = cat(1,x2,x(k));
                y2 = cat(1,y2,y(k));
                
                nx0 = nx(k);
                ny0 = ny(k);
                
                nx2 = cat(1,nx2,nx0);
                ny2 = cat(1,ny2,ny0);
                
                frame2 = cat(1,frame2,frame(k));
                
                xT0 = xT2(iC);
                yT0 = yT2(iC);
                iC = iC+1;
                
                % calculate distance from trail
                dnT0 = sqrt((nx0-xT0).^2+(ny0-yT0).^2);
                shufdnT = cat(1,shufdnT,dnT0);
            end
        end
    end
    
    [initD0,t2s0,spotfound0,dwellT0,quadrant0,initO0,initOb0,spotring0,dtraveled0]=...
        getTime2Spot2(mouse2,sess2,shufdnT,x2,y2,nx2,ny2,frame2,xT2,yT2);
    
    k = ~isnan(t2s0);
    
    sinitD(:,ishuf) = initD0(k);
    st2s(:,ishuf) = t2s0(k);
    sspotfound(:,ishuf) = spotfound0(k);
    sdwellT(:,ishuf) = dwellT0(k);
    squadrant(:,ishuf) = quadrant0(k);
    sinitO(:,ishuf) = initO0(k);
    sinitOb(:,ishuf)=initOb0(k);
    sspotring(:,ishuf)=spotring0(k);
    sdtraveled(:,ishuf)=dtraveled0(k);
    
    fprintf('shuffle %d/%d complete \n',ishuf,nshuf);
    toc;
    
end
