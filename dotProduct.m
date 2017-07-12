for iS=1:3
    k = zeros(length(x0),1); 
    for iL=1:length(xS{iS}); 
        k(round(x0)==xS{iS}(iL) & round(y0)==yS{iS}(iL))=1; 
    end
    k0{iS}=k;
    nT=sum(k);
    F = find(k==1);
    dp{iS}=nan(size(x0));
    vec0 = [vec.x0{iS}./vec.mag{iS},vec.y0{iS}./vec.mag{iS}];
    dp{iS}(F)=dot(repmat(vec0,[length(F),1]),[dx(F),dy(F)],2); 
end