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
d = 0.2;
postSmoothing = 0.1;
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
    notedge = x > 44.8 & x < (1280-44.8) & y > 44.8 & y < (1024-44.8);
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

% JBH comment: added initial velocity maximum threshold below 

V_thres = 500; % 500 cm/s
dx = foaw_diff(x0,dT,m,d,postSmoothing);
dy = foaw_diff(y0,dT,m,d,postSmoothing);
V = sqrt(dx.^2+dy.^2)./11.2; % cm/s
k = V < V_thres & dx < 1000 & dy < 1000;
x0(~k)=nan;
y0(~k)=nan;
nx0(~k)=nan;
ny0(~k)=nan;

% % % %
dnx = foaw_diff(nx0,dT,m,d,postSmoothing);
dny = foaw_diff(ny0,dT,m,d,postSmoothing);
nV = sqrt(dnx.^2+dny.^2)./11.2; % cm/s
k = V < V_thres & dx < 1000 & dy < 1000;
nx0(~k)=nan;
ny0(~k)=nan;


nS = ceil(postSmoothing/dT);
x0(notedge)= medfilt1(x0(notedge),nS,'omitnan','truncate');
y0(notedge)= medfilt1(y0(notedge),nS,'omitnan','truncate');
nx0(notedge)= medfilt1(nx0(notedge),nS,'omitnan','truncate');
ny0(notedge)= medfilt1(ny0(notedge),nS,'omitnan','truncate');

% remove discontinuous 'jumps' in data
k=1;

% JBH comment: reduced maximum jump size to 35 px
discont_thres = 35; % px 35
while any(k)
    % Body jumps
    fx = find(~isnan(x0));
    dfx = diff(x0(fx))./diff(fx);
    dfy = diff(y0(fx))./diff(fx);
    VVb = sqrt(dfx.^2+dfy.^2);
    k = VVb > discont_thres;
    x0(fx(k))=nan;
    y0(fx(k))=nan;
    nx0(fx(k))=nan;
    ny0(fx(k))=nan;
    % JBH comment: performed same comp on nose jumps 
    fx = find(~isnan(nx0));
    dfx = diff(nx0(fx))./diff(fx);
    dfy = diff(ny0(fx))./diff(fx);
    VVn = sqrt(dfx.^2+dfy.^2);
    k = VVn > discont_thres;
    x0(fx(k))=nan;
    y0(fx(k))=nan;
    nx0(fx(k))=nan;
    ny0(fx(k))=nan;
end

% x0(notedge) = cmddenoise(x0(notedge),'db2',2,'s');
% y0(notedge) = cmddenoise(y0(notedge),'db2',2,'s');
% nx0(notedge) = cmddenoise(nx0(notedge),'db2',2,'s');
% ny0(notedge) = cmddenoise(ny0(notedge),'db2',2,'s');

if doTest
    figure(2); clf;
    plot(x,y,'k.'); hold on;
    scatter(x0,y0,[],nanzscore(log10(V)),'filled');
    caxis([-2 2]);
end

end

