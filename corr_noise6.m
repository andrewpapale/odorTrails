function [conc_grid] = corr_noise4(x,y,varargin)
% Add noise, zero pixels for intermittency, then filter
%Inputs: x, y, are vectors of coordinates in the x and y directions. E.g., 
%when the table has dimensions 114x91 cm, x and y take values 
% x = -57:0.1:57;
% y = -45.5:0.1:45.5;

%% Default parameters
sx = 0;             % source x location
sy = 0;             % source y location
stdev = 20;         % spread of exponential
noiselevel = 0.01;  % kc
filt_size = 2;      % Gaussian filter size, in mm
% int_0 = 0.8575;   
% k_int = -0.07691;
int_0 = 1;          % I_0, the intermittency at the source
k_int = -0.02;      % % decay constant of intermittency, in cm^-1

%% Parameter inputs
p = inputParser;
addOptional(p,'sx',sx );
addOptional(p,'sy',sy );
addOptional(p,'sigma',stdev );
addOptional(p,'noise',noiselevel );
addOptional(p,'filt',filt_size );
addOptional(p,'kint',k_int );

parse(p,varargin{:});

sx          = p.Results.sx;
sy          = p.Results.sy;
stdev       = p.Results.sigma;
noiselevel  = p.Results.noise;
filt_size   = p.Results.filt;
k_int       = p.Results.kint;

%% Generate grid

% Define smooth gradient as exponential function
odorFun = @(x,y,varargin)   exp( -sqrt( (x-sx).^2 + (y-sy).^2 ) ./ (2*stdev) );

% Generate grid to evaluate odor function
[X,Y]=meshgrid(x,y);

% Eval function, get smooth gradient
conc_grid = odorFun(X,Y)';

%Generate multiplictative uniform noise of size(conc_grid)
noise_grid = conc_grid.*2.*noiselevel.*(rand(size(conc_grid))-0.5);

% Compute distance from spot on all grid points
d_grid = sqrt((X-sx).^2 + (Y-sy).^2 );

%Compute intermittency value at all grid points based on distance
int_grid = int_0.*exp(k_int.*d_grid');

% Add noise to smooth gradient
conc_grid = conc_grid + noise_grid;

% Set grid points to zero according to intermittency
conc_grid = conc_grid.*(int_grid > rand(size(int_grid)));

% Apply Gaussian filter to resulting grid (note: filter size 0 := no
% smoothing
if filt_size > 0
    conc_grid = imgaussfilt(conc_grid,filt_size);
end


