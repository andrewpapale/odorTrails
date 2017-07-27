function C = Tortuosity1(dx,dy,Ts,m,d,postSmoothing)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



assert(length(dx)==length(dy),'x and y must be the same length');

% [ddx,minMSEx,xSelect] = dxdt(dx,dT,window,postSmoothing);
% [ddy,minMSEy,ySelect] = dxdt(dy,dT,window,postSmoothing);

ddx = foaw_diff(dx, Ts, m, d,postSmoothing);
ddy = foaw_diff(dy, Ts, m, d,postSmoothing);

curv = (dx.*ddy-dy.*ddx);
arcL = (dx.^2+dy.^2).^2;
C = abs(curv)./arcL;




end
