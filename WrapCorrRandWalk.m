% WrapCorrRandWalk

nI = 1;
nM = 5;
nS = 197;
velocity = [1 2 3 4 5 6 7 8 9 10]; % cm/s
nV = length(velocity);
max_turn = [pi/10, pi/9 pi/8 pi/7 pi/6 pi/5 pi/4 pi/3 pi/2 pi]; % rad
nT = length(max_turn);
iC=1;
success = nan(nI*861,nV,nT);
t_final = nan(nI*861,nV,nT);
d_final = nan(nI*861,nV,nT);
d_init = nan(nI*861,nV,nT);
for iM=1:nM
    for iS=1:nS
        kT = mouseT==iM & sessT==iS;
        if sum(kT)>0
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            if ~(xT<44.8 | xT>1280-44.8 | yT<44.8 | yT>1024-44.8) %#ok<OR2>
                x_spot = (xT-(1280/2))/11.2;
                y_spot = (yT-(1040/2))/11.2;
                for iN=1:nI
                    rng('shuffle');
                    for iV=1:nV
                        for iT=1:nT
                            [s0,t0,d0,di] = corrRandWalk(x_spot,y_spot,velocity(iV),max_turn(iT));
                            success(iC,iV,iT)=s0;
                            t_final(iC,iV,iT)=t0;
                            d_final(iC,iV,iT)=d0;
                            d_init(iC,iV,iT)=di;
                        end
                    end
                    fprintf('%d/%d \n',iC,861*nI);
                    iC=iC+1;
                end
                
                
            end
        end
    end
end
