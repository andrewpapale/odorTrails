function [ t_final, pos_vals, nose_vals,dists, th_vals,dirs_vals,success,success_idx,pathlength] = concTrackFun_smooth_comp_ou_p( side_l, conc_thres, capture_dist , varargin)  % vel    vel_min     noiselevel     maxt

%%% Jim Hengenius, August 29, 2016
%%% Mouse odor tracking simulation using comparison of times since last
%%% odor encounter. Consists primarily of Tim Whalen's original tracker
%%% with the addition of scalable arena size, path plotting, and variable
%%% velocity.

%% Options for model variants and behaviors

% Set to one for live animation of mouse position
plotting = 0;

% Set to one for correlated random walk (noise in theta)
randwalk = 1;

% Set to one for odor thresholding, otherwise disregard provided threshold
thresholding = 1;

% Bounded arena? (Reflective walls if 1)
bounded = 1;

if ~thresholding
    conc_thres = 0;
end

% Use real arena dimensions (yes if 1)
arena = 0;
side_lx = 121.920/2;
side_ly = 91.400/2;

% End on capture
success = 0;

% Set to one to linearly decrease body velocity as a function of
% concentration
velocity = 1;
% vel_reduction = 0.5386;

% Capture time t_final
capture = 1;

% Save movie frames?
makemov = 0;


%% Starting postion and heading

% Uncomment following for nearby starting location
% Domain size is 2*side_l x 2*side_l
% side_l = 60;

% X,Y initial conditions

if ~arena
    %     startx = -side_l + 2*rand()*side_l; % Random starting position
    startx = capture_dist + 10 + rand()*(side_l-capture_dist-10);
    starty = -side_l + 2*rand()*side_l;
    theta = 2*rand()*pi; % random starting dirsection
else
    l = 5;
    startx = capture_dist + 2*l + rand()*(side_ly-capture_dist-2*l); % Random starting position
    
    
    
    %Heading initial conditions
    theta = 2*rand()*pi; % random starting dirsection
    starty = 0;
end
%     theta = pi/2;

% If the time measured is for elapsed time unto capture

d = 0;
if capture
    %Reject starting locations w/in cap radius
    d = sqrt(startx^2 + starty^2);
    while d < capture_dist
        startx = -side_l + 2*rand()*side_l; % Random starting position
        % starty = -side_l + 2*rand()*side_l;
        starty = 0;
        d = sqrt(startx^2 + starty^2);
    end
    %Uncomment for theta pointing toward source
    theta = atan2(-starty,-startx );
    thstart = angdiff(theta, atan2(-starty,-startx ));
    
end


%% Mouse intrinsic parameters
% Length in cm, time in sec, angle in radians

vel = 18; % cm/s
vel_min = 0;
l = 5; % radius of perturbation (neck length)
% phi_max = pi/(4); % angle of perturbation NOT SCALED
phi_max = 2*pi/5; % angle of perturbation NOT SCALED
phi = phi_max;
turn = 1; % factor of phi to turn by NOT SCALED
step = .1; % corresponds to 10 Hz sniffing
realtime = 0; % 1 if simulation runs at 1 sec/sec, 0 if runs quickly
internostril_dist = 0.18; % cm
% internostril_dist = 1; % cm
lhat = sqrt(l^2 + (internostril_dist/2)^2);
bin_ang = asin((internostril_dist/2)/lhat);
maxt = 300;% Max sim time
dirs = 0;
noisemod = 1; % Scales noise in heading theta
noiselevel =.01;
phi_effective = 0.3; % Standard deviation of nose deflection
tau = .3; % Time constant for nose deflection
k_bin = 0; % Strength of binaral bias in nose OU random walk

if length(varargin) ==2
    %     wl = varargin{1};
    %     noisemod = varargin{1};
    vel = varargin{1};
    maxt = varargin{2};
    
elseif length(varargin) ==3
    
    vel = varargin{1};
    maxt = varargin{2};
    noiselevel = varargin{3};
elseif length(varargin) ==4
    
    vel = varargin{1};
    vel_min = varargin{2};
    noiselevel = varargin{3};
    maxt = varargin{4};
elseif length(varargin) ==5
    
    vel = varargin{1};
    vel_min = varargin{2};
    noiselevel = varargin{3};
    maxt = varargin{4};
    tau = varargin{5};
elseif length(varargin) ==6
    
    vel = varargin{1};
    vel_min = varargin{2};
    noiselevel = varargin{3};
    maxt = varargin{4};
    tau = varargin{5};
    k_bin = varargin{6};
    
else
    %     wl = 1;
    
end

wl = 1;
wr = 1;
wsum = wl+wr;

% Noise in heading
beta = noisemod*pi/30;
%         beta = .2;
%     beta = rand()*phi_max;
%     noise = beta/phi_max;

%Noise in PID signal
N_pid = 0.10;
%     N_pid = rand();
noise = N_pid;

% Noise in periodic clock
%     gamma = step/100;
gamma = 0;
% Periodic cast amplitude and frequency
amp = phi; % Amplitude of nose arc
freq_max = 1.2; % Frequency of nose cast in seconds^-1
freq = freq_max;


%% Odor Environment

% Point source, Bivariate Gaussian
% stdev = 5;
stdev = 10;
%     odorFun = @(x,y,varargin) odorant_C(0,0,stdev,stdev,0,x,y); % @ origin
odorFun = @(x,y,varargin) odorant_C(0,0,stdev,stdev,0,x,y) ; % @ origin

ofplot = @(x,y) odorant_C(0,0,stdev,stdev,0,x,y); % @ origin;


%% Initialize Plot

if plotting
    figure
    
    % Plot odor function
    scatter(0,0,[],50,'MarkerEdgeColor','g','MarkerFaceColor','g')
    
    if ~arena
        axis([-side_l,side_l,-side_l,side_l])
        axis square
    else % if arena
        axis([-side_lx,side_lx,-side_ly,side_ly])
        %         axis equal
    end
    hold on
end

%Initialize body and nose positions
mx = startx;
my = starty;
nx = startx+l*cos(theta);
ny = starty+l*sin(theta);

if plotting
    N = plot(nx, ny,'v','MarkerFaceColor','red','MarkerSize',6);
    M = plot(mx, my,'o','MarkerFaceColor','black','MarkerSize',6);
    if ~arena
        axis([-side_l,side_l,-side_l,side_l])
        axis square
    else
        axis([-side_lx,side_lx,-side_ly,side_ly])
        %         axis equal
        
    end
    ax = gca;
end

%% Algorithm

% Initialization: Get conc at nose position at the initial position
conc_old = -Inf; % updated in loop before used
conc_new = odorFun(nx, ny,0);




%Scale velocity to step size
v = vel*step;

% Saves variables
pos_vals = nan(maxt/step,2); % x,y position
nose_vals =  nan(maxt/step,2);
th_vals = nan(maxt/step,1);
dirs_vals = nan(maxt/step,1);
% dist_vals = []; % distance from origin
% vel_vals = []; % velocity at each time
phi_vals = nan(maxt/step,1);
th_rel = nan(maxt/step,1);
dists = nan(maxt/step,1);
conc_vals = nan(maxt/step,1);
vel_vals = nan(maxt/step,1);

%Always cast
cast = 1;
F(maxt/step) = struct('cdata',[],'colormap',[]);
del_c = 0;
success_idx = NaN;
final_countdown = -1;

for t = 1:maxt/step
    
    pos_vals(t,:) = [mx my]; %accumulate coordinates
    nose_vals (t,:)= [nx ny];
    th_vals(t) = theta;
    
    
    dists(t) = sqrt(nx^2 + ny^2);
    
    
    if capture
        
        if norm([nx ny]) < capture_dist & success == 0 % Did mouse reach source?
            t_final = t*step;
            success = 1;
            success_idx = t;
            final_countdown = 10;
        end
        if t == maxt/step
            t_final = NaN;              % If time expired, return NaN
        end
        
    end
    
    if t == maxt/step
        t_final = NaN;              % If time expired, return NaN
    end
    
    if success == 1
        final_countdown = final_countdown-1;
    end
    
    if final_countdown == 0
        break                       % If yes, break
    end
    
    
    % Is nose out of bounds?
    if bounded
        if ~arena
            if abs(nx)>side_l || abs(ny)>side_l  %  outside of box?
                %If so, reverse course
                theta = theta + pi; %reverse course
            end
            
        else
            
            % Check if out of bounds for real arena geometry
            if abs(nx)>side_lx || abs(ny) > side_ly  %  outside of box?
                %         If so, reverse course
                theta = theta + pi; %reverse course
            end
            
        end
    end
    
    
    mx = mx + v*cos(theta);
    my = my + v*sin(theta);
    
    
    if cast
        
        
        %Chose a cast dirsection
        
        
        % dirs = dirs + -step*dirs/tau + phi_effective*randn()*sqrt(step/tau); % Standard OU model for nose movement
        dirs = dirs + -step*dirs/tau + phi_effective*randn()*sqrt(step/tau) + k_bin*step*del_c; % Biased OU model with binaral component
        
        dirs_vals(t) = dirs;
        
        
        % Update nose pos
        nx = mx + l*cos(theta+dirs);
        ny = my + l*sin(theta+dirs);
        
        nxl = mx + lhat*cos(theta+dirs+bin_ang);
        nyl = my + lhat*sin(theta+dirs+bin_ang);
        
        nxr = mx + lhat*cos(theta+dirs-bin_ang);
        nyr = my + lhat*sin(theta+dirs-bin_ang);
        conc_vals(t) = conc_new;
        % Sample Concentration
        if (wl*(odorFun(nxl,nyl,conc_new)+max(0, odorFun(nxl,nyl,conc_new)  )   )      +    wr*odorFun(nxr,nyr,conc_new))     /wsum > conc_thres
            conc_old = conc_new;
            
            %             cl = max(0, odorFun(nxl,nyl,conc_old) + (-noiselevel + rand()*2*noiselevel)   ); % Additive uniform noise
            %             cr = max(0, odorFun(nxr,nyr,conc_old) + (-noiselevel + rand()*2*noiselevel)   ); % Additive uniform noise
            
            %             cl = max(0, odorFun(nxl,nyl,conc_old) * (1+noiselevel*randn())   ); %Multiplicative normal noise
            %             cr = max(0, odorFun(nxr,nyr,conc_old) * (1+noiselevel*randn())   ); %Multiplicative normal noise
            
            cl = max(0, odorFun(nxl,nyl,conc_old) * (1+ 2*noiselevel*(randn()-0.5)   )   ); %Multiplicative uniform noise
            cr = max(0, odorFun(nxr,nyr,conc_old) * (1+ 2*noiselevel*(randn()-0.5)   )   ); %Multiplicative uniform noise
            
            conc_new = (wl*cl        +         wr*cr )/wsum;
            
            del_c = (wl*cl - wr*cr)/wsum;
            
            dirs_vals_mc_p(t) = conc_new;
            
            if conc_new>conc_old
                if randwalk
                    theta = theta +dirs*turn+ beta*randn();
                else
                    theta = theta +dirs*turn;
                end
                
            else
                if randwalk
                    theta = theta + beta*randn();
                end
                
            end
            
            
        else % If subthreshold
            conc_old = conc_new;
            conc_new = (wl*odorFun(nxl,nyl,conc_old)+wr*odorFun(nxr,nyr,conc_old))/wsum;
            del_c = (wl*odorFun(nxl,nyl,conc_old) - wr*odorFun(nxr,nyr,conc_old));
            if randwalk
                theta = theta + beta*randn();
            end
        end
        
    end
    phi_vals(t) = dirs;
    
    if double(get(gcf,'CurrentCharacter'))==27
        t_final = t*step;
        break
    end
    
    
    conc = conc_new;
    
    
    % Variable velocities
    if velocity
        hill = 3;
        %         v = (-(vel-vel_min)*conc_new + vel  )*step; % % vel scales from 2 to 10 cm/s (high to low conc)
        v = (vel - (vel-vel_min)*conc_new^hill/(0.5+conc_new^hill))*step;
        
    end
    
    vel_vals(t) = v/step;
    
    
    
    if plotting
        % update graphics
        M.XData = mx;
        M.YData = my;
        N.XData = nx;
        N.YData = ny;
        %         dirs.XData = dx;
        %         dirs.YData = dy;
        
        if arena
            axis([-side_lx,side_lx,-side_ly,side_ly])
            %             axis equal
        end
        
        if realtime
            pause(.1);
        end
        %pause
        drawnow
        if makemov
            F(t) = getframe(ax);
        end
    end
    
end

pos_vals = pos_vals(1:t,:); % x,y position
nose_vals = nose_vals(1:t,:); % nose x,y
th_vals = th_vals(1:t,:); % abs heading angles
th_rel = th_rel(1:t,:); % abs heading angles
dists = dists(1:t,:);


if plotting
    %     for i = 1:t-1
    %
    %         plot(pos_vals(i:i+1,1),pos_vals(i:i+1,2),'k-','LineWidth',1);
    %         plot(nose_vals(i:i+1,1),nose_vals(i:i+1,2),'r-','LineWidth',0.5);
    %     end
    
    
    plot(pos_vals(1:t,1),pos_vals(1:t,2),'k-','LineWidth',1);
    plot(nose_vals(1:t,1),nose_vals(1:t,2),'r-','LineWidth',0.5);
    scatter(nose_vals(1:t,1),nose_vals(1:t,2),50,'r','.');
    
    if arena & bounded
        
        axis([-side_lx,side_lx,-side_ly,side_ly],'equal')
    end
    %
    %     figure;
    %     plot([step:step:t*step],phi_vals(1:t),'.-')
    
    
    if makemov
        F(t) = getframe(ax);
    end
    
    if ~arena
        axis([-side_l,side_l,-side_l,side_l])
        axis square
    elseif bounded
        axis([-side_lx,side_lx,-side_ly,side_ly])
        axis equal
    end
    
end

nose_vals = nose_vals(1:t,:);
pos_vals = pos_vals(1:t,:);
dirs_vals = wrapToPi( dirs_vals(1:t,:) );
th_vals = wrapToPi(   th_vals(1:t,:)   );
th_rel = wrapToPi(   th_rel(1:t,:)   );
thstart = th_vals(1);
if t == maxt/step
    t_final = NaN;              % If time expired, return NaN
else
    t_final = t*step;
end

if success
    pathlength = sum(hypot( diff(nose_vals(:,1)) , diff(nose_vals(:,2)) ) );
else
    pathlength = NaN;
end


t= 0:step:(t-1)/step;

