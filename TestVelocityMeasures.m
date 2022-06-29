% Test Velocity measures

dT = 1/50;
window = 1;
postSmoothing = 0.1;
dN = 0.05;
kp = 10;
ki = 0.1;

% body

dx0 = dxdt(x0,dT,window,postSmoothing);
dy0 = dxdt(y0,dT,window,postSmoothing);
phi = atan2(dy0, dx0);
uphi = unwrap(phi);
dphi0 = dxdt(uphi,dT,window,postSmoothing);

dx1 = foaw_diff(x0,dT,window,dN,postSmoothing);
dy1 = foaw_diff(y0,dT,window,dN,postSmoothing);
phi = atan2(dy1, dx1);
uphi = unwrap(phi);
dphi1 = foaw_diff(uphi,dT,window,dN,postSmoothing);

dx2 = pi_diff(x0,dT,kp,ki,postSmoothing);
dy2 = pi_diff(y0,dT,kp,ki,postSmoothing);
phi = atan2(dy2, dx2);
uphi = unwrap(phi);
dphi2 = pi_diff(uphi,dT,kp,ki,postSmoothing);

dx3 = euler_diff(x0,dT,postSmoothing);
dy3 = euler_diff(y0,dT,postSmoothing);
phi = atan2(dy3,dx3);
uphi = unwrap(phi);
dphi3 = euler_diff(uphi,dT,postSmoothing);

% nose

dnx0 = dxdt(nx0,dT,window,postSmoothing);
dny0 = dxdt(ny0,dT,window,postSmoothing);
phi = atan2(dny0, dnx0);
uphi = unwrap(phi);
ndphi0 = dxdt(uphi,dT,window,postSmoothing);

dnx1 = foaw_diff(nx0,dT,window,dN,postSmoothing);
dny1 = foaw_diff(ny0,dT,window,dN,postSmoothing);
phi = atan2(dny1, dnx1);
uphi = unwrap(phi);
ndphi1 = foaw_diff(uphi,dT,window,dN,postSmoothing);

% dnx2 = pi_diff(nx0,dT,kp,ki,postSmoothing);
% dny2 = pi_diff(ny0,dT,kp,ki,postSmoothing);
% phi = atan2(dny2, dnx2);
% uphi = unwrap(phi);
% ndphi2 = pi_diff(uphi,dT,kp,ki,postSmoothing);
% 
% dnx3 = euler_diff(nx0,dT,postSmoothing);
% dny3 = euler_diff(ny0,dT,postSmoothing);
% phi = atan2(dny3,dnx3);
% uphi = unwrap(phi);
% ndphi3 = euler_diff(uphi,dT,postSmoothing);