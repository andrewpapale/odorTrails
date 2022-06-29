% TimeInRing
% 2017-08-07 AndyP
% get the time in each ring around spots

mint2s = 0;
ktime = 1500; % frames, initial condition
radius = linspace(0,100,100);
nR = length(radius);
nM = 5;
nS = 105;
area = nan(nR,1);
for iR=1:nR
    if iR==1
        area(iR)=pi*radius(iR).^2;
    else
        area(iR)=pi*(radius(iR).^2-radius(iR-1).^2);
    end
end
R = nan(nR,nM,nS);
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            k0 = k & ~edge & mL<25;
            temp = histcounts(dnT(k0),linspace(0,nR,nR));
            R(1:length(temp),iM,iS)=temp;%./nansum(temp,2);
        end
    end
    disp(iM);
end
radii = radius;