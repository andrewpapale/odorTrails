function C = Tortuosity(dx,dy,dT,window,postSmoothing)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



assert(length(dx)==length(dy),'x and y must be the same length');

ddx = dxdt(dx,dT,window,postSmoothing);
ddy = dxdt(dy,dT,window,postSmoothing);

nP = length(dx);

C = nan(nP,1);
for iP=1:nP
     curv = (dx(iP).*ddy(iP)-dy(iP).*ddx(iP))./(dx(iP).^2+dy(iP).^2).^(1.5);
     arcL = sqrt(dx(iP).^2+dy(iP).^2);
     C(iP) = abs(curv)./arcL;
end



end

