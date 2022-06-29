% inWindZone
% 2019-02-25 AndyP

dw = 10*5.174; % cm 'wind zone'

inW = zeros(size(x));
wP = zeros(size(x));
fx2 = fx0;
fy2 = fy0;

vidH = 576;
vidW = 480;

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            fx1 = fx0(iM,iS);
            fy1 = fy0(iM,iS);
            inW0 = inW(k);
            nx0 = nx(k);
            ny0 = ny(k);
            if fy1 > 350 % top
                inW0(nx0 > (fx1-dw) & nx0 < (fx1+dw) & ny0 > 250)=1;
                if fx1 > 300
                    wP0 = ones(sum(k),1);
                else
                    wP0 = ones(sum(k),1)+7;
                end
                fy2(iM,iS)=vidH;
            elseif fx1 > 450 % right
                inW0(ny0 > (fy1-dw) & ny0 < (fy1+dw) & nx0 > 350)=1;
                if fy1 > 250
                    wP0 = ones(sum(k),1)+1;
                else
                    wP0 = ones(sum(k),1)+2;
                end
                fx2(iM,iS)=vidW;
            elseif fy1 < 150 && (fx1 > 100 || fx1 < 450) % bottom
                inW0(nx0 > (fx1-dw) & nx0 < (fx1+dw) & ny0 < 250)=1;
                if fx1 > 300
                    wP0 = ones(sum(k),1)+3;
                else
                    wP0 = ones(sum(k),1)+4;
                end
                fy2(iM,iS)=1;
            elseif fx1 < 100 % left
                inW0(ny0 > (fy1-dw) & ny0 < (fy1+dw) & nx0 < 200)=1;
                if fy1 < 250
                    wP0 = ones(sum(k),1)+5;
                else
                    wP0 = ones(sum(k),1)+6;
                end
                fx2(iM,iS)=1;
            else
            end
            
            inW(k) = inW0;
            wP(k) = wP0;
        end
    end
    disp(iM);
end
            
                
            