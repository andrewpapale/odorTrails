ot = dir('*-Odor-Trail.mat');
pt = dir('*positions.mat');
st = dir('*-StartFrame.mat');

assert(length(ot)==length(pt),'error mismatch in # of files');
assert(length(ot)==length(st),'error mismatch in # of files');

nT = length(ot);
onan = [];
pMM = [];
pmL = [];
pnV = [];
Nmm = [];
for iT=1:nT
    load(ot(iT).name,'xT','yT');
    load(pt(iT).name,'position_results');
    load(st(iT).name,'startFrame');
    
    [x0,y0,nx0,ny0,V,nV0,Nmm0,mL0] = Process_VT(position_results,startFrame);
    nx = position_results.nosePOS(:,1);
    ny = position_results.nosePOS(:,2);
    
    %// Replace NaNs with the closest non-NaNs
    %// Index array for factor
    x1 = 1:numel(x0);
    %// Indices of NaNs
    t2 = find(~isnan(x0));
    nearX = interp1(x1(t2),x0(t2),x1,'nearest');
    nearY = interp1(x1(t2),y0(t2),x1,'nearest');
    edge = ~(nearX>50 & nearY>50 & nearX<nanmax(nearX)-50 & nearY<nanmax(nearY)-50);
    
    % conditionals
    onan0 = sum(isnan(nx(~edge)));
    k = ones(size(nx0));
    k(Nmm0 < 0.5*nanmean(Nmm0) | isnan(Nmm0) & ~edge') = 0;
    pMM0 = sum(k==0);
    k = ones(size(nx0));
    k(mL0>8 & ~edge)=0;
    pmL0 = sum(k==0);
    dnx = foaw_diff(nx0, 1/50, 1, 0.5, 0.2);
    dny = foaw_diff(ny0, 1/50, 1, 0.5, 0.2);
    k = nV0 < 200 & dnx < 1000 & dny < 1000 & ~edge';
    pnV0 = sum(k==0);
    
%     clf; plot(nx,ny,'k.'); hold on;
%     scatter(nx0,ny0,10,nV0,'filled'); colorbar; caxis([0 20]);
%     plot(nx0,ny0,'g.-');
%     plot(yT,xT,'r.','markersize',20);
%     title(strcat(mat2str(iT),'/',mat2str(length(ot))));
     fprintf('%s,%d/%d \n',pt(iT).name,iT,length(ot));
    
    onan = cat(1,onan,onan0);
    pMM = cat(1,pMM,pMM0);
    pmL = cat(1,pmL,pmL0);
    pnV = cat(1,pnV,pnV0);
    Nmm = cat(1,Nmm,0.5*nanmean(Nmm0));
    
    %pause;
end