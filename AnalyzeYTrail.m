% Analyze Y Trail

dtime = 1/50;
m = 50;
d = 1;
postSmoothing = 0.2;

[x0,y0,nx0,ny0] = Process_VT(position_results,startFrame);

xT0 = yT;
yT0 = xT;

xT1 = xT0; %cat(1,xT1,xT0);
yT1 = yT0; %cat(1,yT1,yT0);

%// Index array for factor
x1 = 1:numel(x0);

%// Indices of NaNs
t2 = find(~isnan(x0));

%// Replace NaNs with the closest non-NaNs
nearX0 = interp1(x1(t2),x0(t2),x1,'nearest');
nearY0 = interp1(x1(t2),y0(t2),x1,'nearest');

nearX0(1)=nearX0(2);
nearY0(1)=nearY0(2);

nearX = nearX0';% = cat(1,nearX,nearX0');
nearY = nearY0';% = cat(1,nearY,nearY0');

x = x0; %cat(1,x,x0);
y = y0; %cat(1,y,y0);
nx = nx0; %cat(1,nx,nx0);
ny = ny0; %cat(1,ny,ny0);

frame = (1:length(x0))';

dx = foaw_diff(x0, dtime, m, d, postSmoothing);
dy = foaw_diff(y0, dtime, m, d, postSmoothing);
V = sqrt(dx.^2+dy.^2)./11.2;

ndx = foaw_diff(nx0, dtime, m, d, postSmoothing);
ndy = foaw_diff(ny0, dtime, m, d, postSmoothing);
nV = sqrt(ndx.^2+ndy.^2)./11.2;

% calculate distance from trail
nP = length(x0);
dT = nan(nP,1);
I = nan(nP,1);
dnT = nan(nP,1);
In = nan(nP,1);
for iP=1:nP
    [dT(iP),I(iP)] = nanmin(sqrt((x0(iP)-xT0).^2+(y0(iP)-yT0).^2));
    [dnT(iP),In(iP)] = nanmin(sqrt((nx0(iP)-xT0).^2+(ny0(iP)-yT0).^2));
end
dT = dT./11.2;
Xp = xT0(I);
Yp = yT0(I);
dnT = dnT./11.2;
Xnp = xT0(In);
Ynp = yT0(In);

C = Tortuosity1(dx,dy,dtime,m,d,postSmoothing);
nC = Tortuosity1(ndx,ndy,dtime,m,d,postSmoothing);
idphi = zIdPhi1(dx,dy,dtime,m,d,postSmoothing);
nidphi = zIdPhi1(ndx,ndy,dtime,m,d,postSmoothing);

zidphi= nanzscore(abs(idphi));
znidphi= nanzscore(abs(nidphi));

zC = nanzscore(abs(C));
znC = nanzscore(abs(nC));

dTc = nan(length(x0),1);
dnTc = nan(length(x0),1);
for iP=1:nP
    dTc(iP) = sqrt((x0(iP)-xC).^2+(y0(iP)-yC).^2);
    dnTc(iP) = sqrt((nx0(iP)-xC).^2+(ny0(iP)-yC).^2);
end
dTc = dTc./11.2;
dnTc = dnTc./11.2;

dp = dotProduct(x0,y0,nx0,ny0,xC0,yC0,xC,yC);

dotA = dp{1};
dotB = dp{2};
dotC = dp{3};

inR = inside;

% get angles between trail vectors
vec1 = [vec.x0{1}./vec.mag{1},vec.y0{1}./vec.mag{1}];
vec2 = [vec.x0{2}./vec.mag{2},vec.y0{2}./vec.mag{2}];
vec3 = [vec.x0{3}./vec.mag{3},vec.y0{3}./vec.mag{3}];

thetaA = acos(dot(vec1,vec2))*180/pi;
thetaB = acos(dot(vec2,vec3))*180/pi;
thetaC = acos(dot(vec1,vec3))*180/pi;

theta1 = repmat(thetaA,[length(x0),1]);
theta2 = repmat(thetaB,[length(x0),1]);
theta3 = repmat(thetaC,[length(x0),1]);



temp1 = nan(length(x0),1);
temp2 = nan(length(x0),1);
for iD=1:3
    for iP=1:length(x0)
        
        temp1(iP) = nanmin(sqrt((x0(iP)-xT0(arm{iD}==1)).^2+(y0(iP)-yT0(arm{iD}==1)).^2));
        temp2(iP) = nanmin(sqrt((nx0(iP)-xT0(arm{iD}==1)).^2+(ny0(iP)-yT0(arm{iD}==1)).^2));
        
    end
    dTarm{iD} = temp1./11.2;
    dnTarm{iD} = temp2./11.2;
end

