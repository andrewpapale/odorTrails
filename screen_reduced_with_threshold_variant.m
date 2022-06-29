% Smoothing parameters for dphi calc
m = 20;
d = 0.2;
dtime = 1/10; % Hz
postSmoothing = 0.2; % s


% Compute threshold of 30 cm (300 mm)
side_lx = 114.3/2;
side_ly = 91.400/2;
xgrid = -side_lx:0.1:side_lx;
ygrid = -side_ly:0.1:side_ly;
cgrid = corr_noise6(xgrid,ygrid, 'noise',0.5,'filt',4, 'kint', -0.02, 'sx',-side_lx,'sy',0);
z = cgrid(:,457);
fit = createFit(z);
thres = feval(fit,300);


% Number of parameter sets (50,000) and iterations per set (100)
npars = 50000;
iters = 100;


% Randomize parameter set values
rng('shuffle');
Vmax    = 20  + 20*rand(npars,1);       % [20 40] Velocity max (25)
Vmax(1) = 25;

rng('shuffle');
L       = 0.01+ 9.99*rand(npars,1);     % [0.01 10] Nose length (5)
L(1)    = 5;

rng('shuffle');
dnares  = 0.01+ 0.99*rand(npars,1);     % [0.01 1] Internares distances (0.18)
dnares(1)= 0.18;

rng('shuffle');
sigmin  = 0.01+ 0.99*rand(npars,1);     % [0.01 1] Minimum OU sigma (0.2)
sigmin(1)= 0.2;

rng('shuffle');
sigmax  = 0.01+ 0.99*rand(npars,1);     % [0.01 1] Max OU sigma (0.3)
sigmax(1)= 0.3;

rng('shuffle');
tau  = 0.01+ 0.99*rand(npars,1);     % [0.01 1] OU tau (0.1)
tau(1)= 0.1;

rng('shuffle');
kv    = rand(npars,1);       % [0 1] Concentration k for velocity (0.5)
kv(1) = 0.5; 

rng('shuffle');
kc    = rand(npars,1);    % [0 1] Concentration k for casting (0.5)
kc(1) = 0.5;

rng('shuffle');
nv    = 5*rand(npars,1);    % [0 5] Sigmoid coefficient for velocity (4)
nv(1) = 4;

rng('shuffle');
nc    = 5*rand(npars,1);             % [0 5] Sigmoid coefficient for casting (1)
nc(1) = 1;

rng('shuffle');
ksm    = 5*rand(npars,1);             % [0 5] Binaral magnitude term (1)
ksm(1) = 1;

rng('shuffle');
ksd    = 300*rand(npars,1);             % [0 300] Binaral difference term (200)
ksd(1) = 200;


% Make sure sigmin < sigmax
rng('shuffle');
for i = 1:npars/2
    sigmax(i)  = sigmin(i) + (1-sigmin(i))*rand();  % [sigmin 1]
end

rng('shuffle');
for i = npars/2+1:npars
    sigmin(i)  =0.01+ (sigmax(i)-0.01)*rand();  % [sigmin 1]
end



pars = [Vmax L dnares sigmin sigmax tau kv kc nv nc ksm ksd];


% Vars for storing values
successes       = nan(npars,1);
meanlogdphinV   = nan(npars,1);
meannV          = nan(npars,1);




% Main loop
for i = 1:size(pars,1)
    
    tic
    successlocal = nan(iters,1);
    dphilocal = nan(iters,1);
    nVlocal = nan(iters,1);
    
    parfor j = 1:iters
        
        % Run model iters times
        [~, ~, nose_vals,dists, ~,~,success,~,~,~,~,~] = concTrackFun_smooth_comp_ou_active_nav_variant(100,thres,1.5,...
            'Time', 30,...
            'Ksm',ksm(i),...
            'Ksd',ksd(i),...
            'Vmax', Vmax(i),...
            'L', L(i), ...
            'dnares', dnares(i) ,...
            'phi_min',sigmin(i) ,...
            'phi_max',sigmax(i),...
            'kv', kv(i),...
            'kc', kc(i),...
            'nv', nv(i),...
            'nc', nc(i),...
            'noiselevel', 0.5,...
            'filt',4,...
            'kint', -0.02);
        
        
        % Compute dphi, nV, and log(dphi/nV)
        dphi = abs([0; zIdPhi1(diff(nose_vals(:,2)), diff(nose_vals(:,1)),dtime,m,d,postSmoothing)]);
        nV = [NaN; hypot(diff(  nose_vals(:,1)  ), diff(      nose_vals(:,2)       ))]/dtime;
        logdphinV = log10(dphi./nV);
        
        % Record mean values for dist <= 30 cm
        successlocal(j) = success;
        dphilocal(j)    = nanmean(logdphinV(dists <= 30));
        nVlocal(j)      = nanmean(nV(dists <= 30));
        
        
        
        
    end
    
    % Record mean values across iterations
    successes(i)        = nanmean(successlocal);
    meanlogdphinV(i)    = nanmean(dphilocal);
    meannV(i)           = nanmean(nVlocal);
    
    if mod(i,100)==0
        disp([num2str(i) ' of ' num2str(length(pars)) '.  Fraction complete: ' num2str(i/length(pars)) '. Success rate: ' num2str(successes(i))])
        toc
    end
    
end

% Clear temp vars
clear dphi nV logdphinV successlocal dphilocal nVlocal m d  dtime  postSmoothing  cgrid z fit thres i j side_lx side_ly xgrid ygrid
    