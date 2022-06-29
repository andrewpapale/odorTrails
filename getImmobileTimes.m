function [sumTimmobile,immobile,V,nV] = getImmobileTimes(x0,y0,nx0,ny0,arena_data,thresholdbody,thresholdnose,bins,timethreshold)
% 2017-10-13 AndyP
% 2017-10-30 AndyP, modified to use second threshold based on nose position
% get times that the animal is immobile
% [sumTimmobile,immobile,V,nV] = getImmobileTimes(x0,y0,nx0,ny0,arena_data,thresholdbody,thresholdnose,,bins,timethreshold)
% Example: [sumTimmobile,immobile,V] = getImmobileTimes(x0,y0,nx0,ny0,arena_data,1,0.5,[1 300 600],60);

doTest = true;

m = 30;
d = 0.5;
postSmoothing = 1;
dT = 1/30;
conv_factor = arena_data.pixels_per_mm*10;

dx = foaw_diff(x0,dT,m,d,postSmoothing);
dy = foaw_diff(y0,dT,m,d,postSmoothing);
V = sqrt(dx.^2+dy.^2)./conv_factor; % cm/s

dnx = foaw_diff(nx0,dT,m,d,postSmoothing);
dny = foaw_diff(ny0,dT,m,d,postSmoothing);
nV = sqrt(dnx.^2+dny.^2)./conv_factor; % cm/s

immobile = nan(size(V));
for iT=timethreshold:timethreshold:length(V)    
    meanVpertime = nanmean(V(iT-timethreshold+1:iT));
    meannVpertime = nanmean(nV(iT-timethreshold+1:iT));
    if meanVpertime < thresholdbody & meannVpertime < thresholdnose %#ok<AND2>
        immobile(iT-timethreshold+1:iT)=1;
    else
        immobile(iT-timethreshold+1:iT)=0;
    end
end

sumTimmobile = histcn(x0,bins,'AccumData',immobile,'fun',@nansum);
sumTimmobile = sumTimmobile/30;

if doTest
   scatter(x0,y0,20,V);
   hold on;
   plot(x0(immobile==1),y0(immobile==1),'rx');
end

end

