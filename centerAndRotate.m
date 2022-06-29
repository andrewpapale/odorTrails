% centerAndRotate
% center and rotate
% 2019-03-01 AndyP

nxc = nx-cx1;
nyc = ny-cy1;

xc = x-cx1;
yc = y-cy1;

fx2 = nan(size(x));
fy2 = nan(size(x));

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            k1 = k & ~isnan(nx);
            fxtemp = fx0(iM,iS);
            fytemp = fy0(iM,iS);
            switch sP(iM,iS)
                case 1
                    ix = interp1(unique(nx(k1)),unique(nx(k1)),fxtemp,'nearest','extrap');
                    ny0 = ny(k);
                    dS = nanmin(ny0(ix));
                    fy2(k) = repmat(dS,[sum(k),1]);
                    fx2(k) = repmat(fxtemp,[sum(k),1]);
                case 2
                    ix = interp1(unique(ny(k1)),unique(ny(k1)),fytemp,'nearest','extrap');
                    nx0 = nx(k);
                    dS = nanmax(nx0(ix));
                    fx2(k) = repmat(dS,[sum(k),1]);
                    fy2(k) = repmat(fytemp,[sum(k),1]);
                case 3
                    ix = interp1(unique(nx(k1)),unique(nx(k1)),fxtemp,'nearest','extrap');
                    ny0 = ny(k);
                    dS = nanmax(ny0(ix));
                    fy2(k) = repmat(dS,[sum(k),1]);
                    fx2(k) = repmat(fxtemp,[sum(k),1]);
                case 4
                    ix = interp1(unique(ny(k1)),unique(ny(k1)),fytemp,'nearest','extrap');
                    nx0 = nx(k);
                    dS = nanmin(nx0(ix));
                    fx2(k) = repmat(dS,[sum(k),1]);
                    fy2(k) = repmat(fytemp,[sum(k),1]);
            end
        end
    end
end



nxf = nx-fx2;
nyf = ny-fy2;

cxc = cx1-cx1;
cyc = cy1-cy1;

cxfs = cx1-fx2;
cyfs = cy1-fy2;

alpha = nan(size(x));
alpha(sP0==1) = 270;
alpha(sP0==2) = 180;
alpha(sP0==3) = 90;
alpha(sP0==4) = 0;

torient = wrapTo360(torient-alpha);

M = [nxc, nyc]';
bM = [xc,yc]';
Mf = [nxf,nyf]';
Ms = [cxc,cyc]';
Fs = [cxfs,cyfs]';
M1 = nan(size(M));
Mf1 = nan(size(M));
M2 = nan(size(M));
M3 = nan(size(M));
for iT=1:length(alpha)
    R = [cosd(alpha(iT)), -sind(alpha(iT)); sind(alpha(iT)), cosd(alpha(iT))];
    M1(:,iT) = R*M(:,iT);
    Mf1(:,iT) = R*Mf(:,iT);
    M2(:,iT) = R*Ms(:,iT);
    M3(:,iT) = R*Fs(:,iT);
end
nxcr = M1(1,:);
nycr = M1(2,:);

tro = [];

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
             xc0 = xc(k);
             yc0 = yc(k);
             nxc0 = nxc(k);
             nyc0 = nyc(k);
             
             BSx = xc0;
             BSy = yc0;
             BNx = nxc0-xc0;
             BNy = nyc0-yc0;
             
             k0 = ~isnan(BSx) & ~isnan(BNx);
             Ik = find(k0==1);
             firstk = find(k==1,1,'first');

%              A = atan2(yT-y0,xT-x0);
%              B = atan2(ny0-y0,nx0-x0);
%             orient(k)=angdiff(A,B)*180/pi;
            torient0 = nan(sum(k),1);
            torient0(k0) = atan2d(BNy(k0),BNx(k0));
%             A = sqrt(BSx.^2+BSy.^2);
%             B = sqrt(BNx.^2+BNy.^2);
%             torient0(k0) = acosd(dot([BSx(k0),BSy(k0)]',[BNx(k0),BNy(k0)]')'./(A(k0).*B(k0)));
             tro = cat(1,tro,torient0);


        end
    end
end

nxfr = Mf1(1,:);
nyfr = Mf1(2,:);

cxcr = M2(1,:);
cycr = M2(2,:);

fxcr = M3(1,:);
fycr = M3(2,:);