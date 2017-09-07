function [x0,y0,nx0,ny0,V,nV,Nmm,mL] = Process_VT1(position_results,startFrame)
% 2017-06-08 AndyP
% 2017-06-09 AndyP, changed condition to vec < P on line 54/line 58, added
% maxIter of 1000 for gmm fit
% 2017-06-15 AndyP, removed gmm, now using simple thresholding and velocity
% thresholding
% [x0,y0,nx0,ny0,V,nV] = Process_VT(position_results);
% processes position coordinates from optimouse to clean up trajectories
%
% INPUTS
% position_results is the position output structure from optimouse
% OUTPUTS
% x0,y0 outputs are the processed body positions
% nx0,ny0 are the processed nose positions
% V,nV is the velocity of the body and nose, respectively

doTest = false;
window = 1;
postSmoothing = 0.1;
dT = 1/50;


x = position_results.mouseCOM(startFrame:end,1);
y = position_results.mouseCOM(startFrame:end,2);

nx = position_results.nosePOS(startFrame:end,1);
ny = position_results.nosePOS(startFrame:end,2);

x0 = x;
y0 = y;

%notedge = x > 50+0.50 & x < (1280.50-50) & y > 50+0.50 & y < (1018.50-50);

nx0 = nx;
ny0 = ny;
% nx0(~notedge)=nan;
% ny0(~notedge)=nan;

mm = position_results.MouseMean(startFrame:end);
k = ones(length(x),1);
mm(mm==0)=nan;
Nmm = log10(mm)./nanmax(log10(mm));
k(Nmm < 0.5*nanmean(Nmm) | isnan(Nmm)) = 0;
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

mL = position_results.mouse_length(startFrame:end)./14.8;
k = ones(size(x0));
k(mL>8)=0; % mice are <8cm long
% x0(~k)=nan;
% y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

nS = ceil(postSmoothing/dT);
x0= medfilt1(x0,nS,'omitnan','truncate');
y0= medfilt1(y0,nS,'omitnan','truncate');
nx0= medfilt1(nx0,nS,'omitnan','truncate');
ny0= medfilt1(ny0,nS,'omitnan','truncate');

dx = dxdt(x0,dT,window,postSmoothing);
dy = dxdt(y0,dT,window,postSmoothing);
V = sqrt(dx.^2+dy.^2)./14.8; % cm/s
k = V < 200 & dx < 1000 & dy < 1000;
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

% % % %
dnx = dxdt(nx0,dT,window,postSmoothing);
dny = dxdt(ny0,dT,window,postSmoothing);
nV = sqrt(dnx.^2+dny.^2)./14.8; % cm/s
k = nV < 200 & dnx < 1000 & dny < 1000;
nx0(~k)=nan;
ny0(~k)=nan;


if doTest
    figure(2); clf;
    plot(x,y,'k.'); hold on;
    scatter(x0,y0,[],nanzscore(log10(V)),'filled');
    caxis([-2 2]);
end

end

