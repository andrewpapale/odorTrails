spotfound1=nan(size(x)); 
spotfound0trunc = nan(size(x));
fpast = 0;
lastt2s = 1500; % frames, initial value
for iM=1:nM
    for iS=1:nS
        if t2s(iM,iS)>0 
            k = mouse==iM & sess==iS; 
            spotfound1(k)=1; 
            lastt2s = round(t2s(iM,iS)*50);
            spotfound0trunc(k)=0;
        elseif t2s(iM,iS)==-20
            k = mouse==iM & sess==iS; 
            spotfound1(k)=0;
            kT = find(k);
            kT1 = kT(kT<=min(lastt2s+kT(1)+fpast,max(kT)));
            kT2 = kT(kT>min(lastt2s+kT(1)+fpast,max(kT)));
            spotfound0trunc(kT1)=1;
            spotfound0trunc(kT2)=0;
        end
    end
    disp(iM); 
end