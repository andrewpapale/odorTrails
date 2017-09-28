function D = ComputeDistanceTraveled(x0,y0)
% 2017-09-27 AndyP
% D = ComputeDistanceTraveled(x0,y0);
% compute distance traveled using body center of mass positions <x0,y0>.

nT = length(x0);
D = nan(nT,1);
for iT=2:nT
    D(iT)=sqrt((x0(iT)-x0(iT-1)).^2+(y0(iT)-y0(iT-1)).^2);
end


end

