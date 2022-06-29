function [success,t_final,d_final,d_init] = corrRandWalk(x_spot,y_spot,velocity,maxTurn)

%Random walk parameters
%velocity = 10; % cm/s
%maxTurn = pi/3; % Diffusion term

doplot = true;

max_time = 300; % seconds
dt = 0.02; % seconds
success = 0; % Flag for successful capture
t_final = NaN;
d_final = NaN;

%Arena parameters
xmin = -60;
xmax =  60;
ymin = -45;
ymax =  45;

%x_spot = 0;
%y_spot = 0;

cap_radius = 2;


% Initial conditions
x_start = xmin + (xmax-xmin)*rand();
y_start = ymin + (ymax-ymin)*rand();
theta_start = -pi + 2*pi*rand();

% Storage arrays and initial conditions

x = nan(max_time/dt,1);
    x(1) = x_start;
y = nan(max_time/dt,1);
    y(1) = y_start;
theta = nan(max_time/dt,1); % Heading
    theta(1) = theta_start;
d_spot = nan(max_time/dt,1);
    d_spot(1) = sqrt((x_start-x_spot).^2 + (y_start-y_spot).^2    );

    d_init = d_spot(1);

%Main loop
cap_ct = 0;
for t = 2:max_time/dt
   
    
    %Update body position
    x(t) = x(t-1) + velocity*dt*cos(theta(t-1));
    y(t) = y(t-1) + velocity*dt*sin(theta(t-1));
    
    % Compute distance to source
    d_spot(t) = sqrt( (x(t) - x_spot).^2 +  (y(t) - y_spot).^2  );
    if d_spot(t) < cap_radius % If successful, break loop
        cap_ct = cap_ct+1;
    else
        cap_ct = 0;
    end
    if cap_ct > 0.5/dt
        success = 1;
        t_final = t*dt;
        d_final = nansum(sqrt(diff(x(1:t)).^2+diff(y(1:t)).^2));
        break
    end
    if x(t) < xmin | x(t) > xmax | y(t) < ymin | y(t) > ymax %#ok<OR2> % Check if out of bounds; if so, reverse course
        theta(t) = theta(t-1)-pi;
%         x(t) = x(t-1);
%         y(t) = y(t-1);
        
    else % Else, update heading
       
            theta(t) = theta(t-1) + sqrt(dt)*randn()*maxTurn;

    end
        
    
end


if doplot
   keyboard; 
end
