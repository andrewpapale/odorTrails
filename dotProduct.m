function dp = dotProduct(x0,y0,nx0,ny0,xC0,yC0,xC,yC)
  
dx = (nx0-x0);
dy = (ny0-y0);
mag = sqrt(dx.^2+dy.^2);
vec1 = [dx./mag,dy./mag];

for iS=1:3
    dp{iS}=nan(size(x0));
    mag = sqrt((xC0{iS}-xC).^2+(yC0{iS}-yC).^2);
    vec0 = [(xC0{iS}-xC)./mag,(yC0{iS}-yC)./mag];
    dp{iS}=dot(repmat(vec0,[length(x0),1]),vec1,2); 
end