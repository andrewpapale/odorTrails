% fixOrientation
temp = orient;
nM = max(mouse);
nS = max(sess);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
             temp(k)= nanfastsmooth(temp(k),5,1,0);
%             mO = nanmedian(temp(k));
%             if mO > -180 && mO < 180 % ok
                 N = 0;
%             else
%                 N = round(mO/180);
%             end
           % temp(k)=wrapTo180(temp(k)-N*180);
        end
    end
    disp(iM);
end

