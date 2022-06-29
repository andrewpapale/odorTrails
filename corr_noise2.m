function [conc_grid] = corr_noise2(x,y,varargin)

sx = 0;
sy = 0;
stdev = 20;
noiselevel = 0.01;
filt_size = 2;

p = inputParser;
addOptional(p,'sx',sx );
addOptional(p,'sy',sy );
addOptional(p,'sigma',stdev );
addOptional(p,'noise',noiselevel );
addOptional(p,'filt',filt_size );


parse(p,varargin{:});

sx          = p.Results.sx;
sy          = p.Results.sy;
stdev       = p.Results.sigma;
noiselevel  = p.Results.noise;
filt_size   = p.Results.filt;

odorFun = @(x,y,varargin)   exp( -sqrt( (x-sx).^2 + (y-sy).^2 ) ./ (2*stdev) ); % exponential


conc_grid = nan(length(x),length(y));
noise_grid = nan(length(x),length(y));

[X,Y]=meshgrid(x,y);

conc_grid = odorFun(X,Y)';
noise_grid = conc_grid.*2.*noiselevel.*(rand(size(conc_grid))-0.5);
 
if filt_size > 0
    noise_grid = imgaussfilt(noise_grid,filt_size);
end

conc_grid = conc_grid + noise_grid;

conc_grid(conc_grid<0) = 0;