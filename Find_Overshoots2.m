% 2018-07-23 AndyP
% Find_Overshoots2

%outTthr = 1; % s
%outTthr = 5; % cm
outTthr1 = linspace(2,100,20);
nO = length(outTthr1);
nM = max(mouse);
nS = max(sess);
nOver = nan(nM,nS,nO);

for iO=1:nO
    outTthr = outTthr1(iO);
    for iM=1:nM
        for iS=1:nS
            if t2s(iM,iS)> 2
                % find entry point
                k = mouse==iM & sess==iS;
                x0 = x(k);
                y0 = y(k);
                nx0 = nx(k);
                ny0 = ny(k);
                %dx = foaw_diff(x0, 1/50, 50, 0.5,0.1);
                %dy = foaw_diff(y0, 1/50, 50, 0.5,0.1);
                dnT0 = dnT(k);
                v0 = 1:sum(k);
                frame0 = round(t2s(iM,iS)*50);
                %outerR = find(v0 > frame0-outTthr*50,1,'first');
                outerR = find(dnT0 < outTthr & v0 < frame0,1,'first');
                
                %dy0 = (y0(outerR)+dy(outerR))-y0(outerR);
                %dx0 = ((x0(outerR)+dx(outerR))-x0(outerR));
                %m0 = dy0./dx0;
                %m0 = -1./m0;
                
                kT = mouseT==iM & sessT==iS;
                %cx0 = nanmedian(xT1(kT));
                %cy0 = nanmedian(yT1(kT));
                
                %result = computeline([cx0,cy0],m0, [1,1300]);
                
                %nP = length(result);
                %x1 = result{1}(:,1);
                %y1 = result{1}(:,2);
                %x2 = result{nP}(:,1);
                %y2 = result{nP}(:,2);
                
                v2 = outerR:frame0;
                
                if any(dnT0(v2)>outTthr)
                    nOver(iM,iS,iO)=1;
                else
                    nOver(iM,iS,iO)=0;
                end
                
                % find on which side the entry point is relative to the
                % bisecting line
                %Dline = (x0(v2)-x1).*(y2-y1)-(y0(v2)-y1).*(x2-x1);
                %Dline = Dline > 0;
                
                %             if any(Dline)~=Dline(1) % then mouse has crossed the line
                %                 nOver(iM,iS)=1;
                %                 %keyboard;
                %             else
                %                 nOver(iM,iS)=0;
                %             end
                
            end
        end
        disp(iM);
    end
    disp(iO);
end