% 2017-07-12 AndyP
% Analyze one session 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

C = [];
nC = [];
idphi = [];
nidphi = [];
mouse = [];
trial = [];
sess = [];
conc = [];
frame = [];
mouseT = [];
sessT = [];


dtime = 1/50;
window = 1;
postSmoothing = 0.5;

xT0 = yT;
yT0 = xT;

xT1 = xT0; %cat(1,xT1,xT0);
yT1 = yT0; %cat(1,yT1,yT0);



[x0,y0,nx0,ny0,~,~] = Process_VT(position_results,startFrame);

% for each nan, find nearest position to each nan...

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

frame = cat(1,frame,(1:length(x0))');

dx = dxdt(x0,dtime,window,postSmoothing);
dy = dxdt(y0,dtime,window,postSmoothing);
V0 = sqrt(dx.^2+dy.^2)./11.2;
V = V0;% cat(1,V,V0); % cm/s

ndx = dxdt(nx0,dtime,window,postSmoothing);
ndy = dxdt(ny0,dtime,window,postSmoothing);
nV0 = sqrt(ndx.^2+ndy.^2)./11.2;
nV = nV0; % cat(1,nV,nV0); % cm/s
% calculate distance from trail
nP = length(x0);
dT0 = nan(nP,1);
I = nan(nP,1);
for iP=1:nP
    [dT0(iP),I(iP)] = nanmin(sqrt((x0(iP)-xT0).^2+(y0(iP)-yT0).^2));
    [dnT0(iP),In(iP)] = nanmin(sqrt((nx0(iP)-xT0).^2+(ny0(iP)-yT0).^2));
end
dT = dT0;
Xp = xT0(I);
Yp = yT0(I);
dnT = dnT0;
Xnp = xT0(In);
Ynp = yT0(In);

C0 = Tortuosity(dx,dy,dtime,window,postSmoothing);
C = cat(1,C,C0);
%
nC0 = Tortuosity(ndx,ndy,dtime,window,postSmoothing);
nC = cat(1,nC,nC0);

idphi0 = zIdPhi(dx,dy,dtime,window,postSmoothing);
nidphi0 = zIdPhi(ndx,ndy,dtime,window,postSmoothing);

idphi = cat(1,idphi,idphi0);
nidphi = cat(1,nidphi,nidphi0);

zidphi=nanzscore(abs(idphi));
znidphi = nanzscore(abs(nidphi));

zC = nanzscore(abs(C));
znC = nanzscore(abs(nC));


