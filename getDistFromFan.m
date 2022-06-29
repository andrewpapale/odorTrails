% getDistFromFan
% 2019-02-27 AndyP

dnF = nan(size(x));

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            dnF(k) = sqrt((nx(k)-fx2(k)).^2+(ny(k)-fy2(k)).^2)./5.174;
        end
    end
    disp(iM);
end