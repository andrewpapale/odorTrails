function [x,y,th,dnT,idphi,v,V]=randwalk_01(frame,vdist,T)

% Discrete-time continuous-space random walk (constant velocity)
% Jim Hengenius - 08/21/2017

xT = T(1);
yT = T(2);

t_max = length(frame)-1; % Elapsed time (frames)
%dt = 0.02; % Time step to increment simulation

vdist = vdist(~isnan(vdist));
indx = randi(length(vdist),[t_max,1]);
vel = vdist(indx); % Velocity (pixels/s)
v = vel; % Scaled for dt

% Boundaries for rectangular domain
y_m = 1020/2;
x_m = 1280/2;

% Initial conditions
x_0 = 1020/2;
y_0 = 1280/2;
th_0 = 3*pi/2; % 270 deg

deltath = pi/4; % Maximum change in heading theta per unit time
dth = deltath;

% Arrays to store values
x = [ x_0 ;  nan(t_max,1) ];
y = [ y_0 ;  nan(t_max,1) ];
th = [ th_0 ;  nan(t_max,1) ]; % Note: t-th timestep uses (t-1)-th th value. Just an convention.


% Main  loop
for t = 2:t_max
    
    % Check if last step took mouse out of bounds (OOB)
    if (abs(y(t-1)) > y_m ||   abs(x(t-1)) > x_m) && t>2
        % If OOB, revert back to last in-bounds position, flip heading 180*
        x(t-1) = x(t-2);
        y(t-1) = y(t-2);
        th(t-1) = th(t-1) + pi;
    end  

    
    % Add noise to heading (applied in next step) - uncomment to choose
    % type of distribution:
    
        % Uniform random angle within dth (correlated random walk):
        th(t) = th(t-1) - dth + 2*dth*rand();
        
        % Normal random angle scaled by dth  (correlated random walk):
%         th(t) = th(t-1) + dth*randn()/6; %Falls w/in dth ~99% of time

        % Uniform random angle on unit circle (uncorrelated random walk)
%         th(t) = -2*pi + 4*pi*rand();
        


    % Update position
    x(t) = x(t-1) + v(t-1).*cos(th(t-1));
    y(t) = y(t-1) + v(t-1).*sin(th(t-1));

end

x = x+1280/2;
y = y+1020/2;

% figure % Plot 2D histogram
% histogram2(x,y,[100 100],'facecolor','flat')
% colorbar
% title('Position Density')
% xlabel('X (cm)')
% ylabel('Y (cm)')
% 
% figure % Rose-plot of headings
% rose(wrapTo2Pi(th))
% title('Heading \theta Distribution')

dx = foaw_diff(x, 1/50, 50, 0.5, 0.2);
dy = foaw_diff(y, 1/50, 50, 0.5, 0.2);
idphi = zIdPhi1(dx,dy,1/50,50,0.5,0.2);
V = sqrt(dx.^2+dy.^2);
% log10idp = logVar(abs(idphi)./V);
% znP = nanzscore(log10idp);

dnT = sqrt((x-xT).^2+(y-yT).^2);

end





