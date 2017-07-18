function dp = dotProduct(x0,y0,nx0,ny0,vec)
  
dx = (nx0-x0);
dy = (ny0-y0);
mag = sqrt(dx.^2+dy.^2);
vec1 = [dx./mag,dy./mag];

for iS=1:3
    dp{iS}=nan(size(x0));
    vec0 = [vec.x0{iS}./vec.mag{iS},vec.y0{iS}./vec.mag{iS}];
    dp{iS}=dot(repmat(vec0,[length(x0),1]),vec1,2); 
end