function [C,minMSEx,xSelect,minMSEy,ySelect] = Tortuosity(dx,dy,dT,window,postSmoothing)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



assert(length(dx)==length(dy),'x and y must be the same length');

[ddx,minMSEx,xSelect] = dxdt(dx,dT,window,postSmoothing);
[ddy,minMSEy,ySelect] = dxdt(dy,dT,window,postSmoothing);

curv = (dx.*ddy-dy.*ddx);
arcL = (dx.^2+dy.^2).^2;
C = abs(curv)./arcL;




end

