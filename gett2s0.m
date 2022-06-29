% get t2s0
t2s0 = zeros(size(x)); 
t2s20 = zeros(size(x)); 
fpast = 0;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS; 
        if sum(k)>0
            if t2s(iM,iS)>1 
                temp = round(t2s(iM,iS)*50);
                first = find(k==1,1,'first'); 
                t2s0(first:temp+first+fpast)=1; 
%                 if t2s2(iM,iS)>1
%                     tstart = round(max(t2s2(iM,iS)*50-500,t2s(iM,iS)*50+fpast));
%                     last = find(k==1,1,'last');
%                     temp = zeros(sum(k),1);
%                     temp(1:min(last,tstart+fpast))=1;
%                     t2s20(k)=temp;
%                 end
            end 
        end 
    end 
    disp(iM);
end