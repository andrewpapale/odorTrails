function [dp,k0] = dotProduct(x0,y0,xS,yS,vec)

dx = dxdt(x0,1/50,1,0.1);
dy = dxdt(y0,1/50,1,0.1);    
rx0 = round(x0);
ry0 = round(y0);
for iS=1:3
    k = zeros(length(x0),1);

    for iL=1:length(xS{iS}); 
        k(rx0==xS{iS}(iL) & ry0==yS{iS}(iL))=1; 
    end
    k0{iS}=k;
    F = find(k==1);
    dp{iS}=nan(size(x0));
    vec0 = [vec.x0{iS}./vec.mag{iS},vec.y0{iS}./vec.mag{iS}];
    dp{iS}(F)=dot(repmat(vec0,[length(F),1]),[dx(F),dy(F)],2); 
end