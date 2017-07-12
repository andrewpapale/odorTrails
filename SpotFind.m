
minD = nan(3,45);
frameM = nan(3,45);
xnM = nan(3,45);
ynM = nan(3,45);
timeInSpot = nan(3,45);

locs = find(dnT./11.2<0.5);

kfound = zeros(size(dnT));
for iM=1:3
    for iS=1:45
        k = find(mouse==iM & sess==iS);
        locs0 = [];
        for iT=1:length(locs)
            if any(k==locs(iT));
                locs0 = cat(1,locs0,locs(iT));
            end
        end
        firstSpot = min(locs0);
        if ~isempty(firstSpot)
            kfound(k)=1;
            minD(iM,iS)=dnT(firstSpot);
            frameM(iM,iS) = frame(firstSpot);
            xnM(iM,iS) = nx(firstSpot);
            ynM(iM,iS) = ny(firstSpot);
            %timeInSpot(iM,iS) = w(find(locs==firstSpot)); unreliable
        end
    end
end

bins = linspace(0,1600,1000);
clear H H1
nanNot = 0;
nanFound  = 0;
Hnot = nan(99,99,3,45);
Hfound = nan(99,99,3,45);
for iM=1:3
    for iS=1:45
        k = find(mouse==iM & sess==iS);
        if ~isempty(k)
            if isnan(minD(iM,iS))
                nanNot = nanNot+sum(isnan(x(k)));
                H0 = hist(dnT(k),bins);
                H(iM,:,iS)=H0;
                H0 = histcn([x(k),y(k)],linspace(1,1300,100),linspace(1,1000,100));
                Hnot(1:size(H0,1),1:size(H0,2),iM,iS)=H0;
            else
                nanFound = nanFound+sum(isnan(x(k)));
                H0 = hist(dnT(k),bins);
                H1(iM,:,iS)=H0;
                H0 = histcn([x(k),y(k)],linspace(1,1300,100),linspace(1,1000,100));
                Hfound(1:size(H0,1),1:size(H0,2),iM,iS)=H0;
            end
        end
    end
end

