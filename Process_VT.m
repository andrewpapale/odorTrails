function [x0,y0,nx0,ny0,V,nV,notedge,Nmm,mL] = Process_VT(position_results,startFrame,edgeflag)
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
m = 50;
d = 0.5;
postSmoothing = 0.2;
dT = 1/50;


temp = position_results.mouseCOM;
x = temp(startFrame:end,1);
y = temp(startFrame:end,2);

temp = position_results.nosePOS;
nx = temp(startFrame:end,1);
ny = temp(startFrame:end,2);

x0 = x;
y0 = y;

if ~edgeflag
    notedge = x > 50 & x < (1280-50) & y > 50 & y < (1024-50);
else
    warning('spot is at the edge for this trial...using truncated border of 1 pixel');
    notedge = x > 1 & x <(1280-1) & y > 1 & y < (1024-1);
end
    
nx0 = nx;
ny0 = ny;
nx0(~notedge)=nan;
ny0(~notedge)=nan;

mm = position_results.MouseMean;
mm = mm(startFrame:end);
k = ones(size(x));
mm(mm==0)=nan;
Nmm = log10(mm)./nanmax(log10(mm));
k(Nmm < 0.5*nanmean(Nmm) | isnan(Nmm)) = 0;
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

mL = position_results.mouse_length;
mL = mL(startFrame:end)./11.2;
k = ones(size(x0));
k(mL>8)=0; % mice are <8cm long
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

nS = ceil(postSmoothing/dT);
x0(notedge)= medfilt1(x0(notedge),nS,'omitnan','truncate');
y0(notedge)= medfilt1(y0(notedge),nS,'omitnan','truncate');
nx0(notedge)= medfilt1(nx0(notedge),nS,'omitnan','truncate');
ny0(notedge)= medfilt1(ny0(notedge),nS,'omitnan','truncate');

dx = foaw_diff(x0,dT,m,d,postSmoothing);
dy = foaw_diff(y0,dT,m,d,postSmoothing);
V = sqrt(dx.^2+dy.^2)./11.2; % cm/s
k = V < 200 & dx < 200 & dy < 200;
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

% % % %
dnx = foaw_diff(nx0,dT,m,d,postSmoothing);
dny = foaw_diff(ny0,dT,m,d,postSmoothing);
nV = sqrt(dnx.^2+dny.^2)./11.2; % cm/s
k = nV < 200 & dnx < 200 & dny < 200;
nx0(~k)=nan;
ny0(~k)=nan;


if doTest
    figure(2); clf;
    plot(x,y,'k.'); hold on;
    scatter(x0,y0,[],nanzscore(log10(V)),'filled');
    caxis([-2 2]);
end

end

