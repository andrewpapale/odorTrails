angcast1 = nan(size(x)); 
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS; 
        if sum(k)>0
            temp=abs(orient(k))>45; 
            k0 = ~edge(k) & ~isstopped(k); 
            dnT0 = dnT(k); 
            [H0,~,~,Hbin]=histcn(dnT0(k0),linspace(0,120,120),'AccumData',temp(k0),'fun',@nansum); 
            k1 = Hbin~=0;
            angcast1(k1)=temp(k1)./H0(Hbin(k1)); 
        end 
    end 
    disp(iM); 
end