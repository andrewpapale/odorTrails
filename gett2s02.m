% get t2s0
t2s0 = zeros(size(x)); 
t2s20 = zeros(size(x));
fpast = 0;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS; 
        kS = mouseT==iM & sessT==iS; % nTrials counter
        if sum(k)>0
            if t2s(kS)>1 
                f2s0 = f2s(kS);
                f2f0 = f2f(kS);
                first = find(k==1,1,'first'); 
                t2s0(first:first+f2s0(1)+fpast)=1; 
                if length(f2s0)>1 && ~(f2s0(2)-f2f0(1)>=f2s0(1)+fpast)
                    t2s20(first+f2f0(1):first+f2s0(2)+fpast)=1;
                end
                if any(t2s20+t2s0==2)
                    keyboard; % debug
                end
            end 
        end 
    end
    disp(iM);
end