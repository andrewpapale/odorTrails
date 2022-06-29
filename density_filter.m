function [xd , yd , rho ] = density_filter(x , y , n_dist, n_thres)
% x and y are coordinates for a nose or body trajectory in cm

% n_dist is the distance threshold used to compute the local density at
% each point (number of neighbors within n_dist of a point)

% n_thres is the threshold number of neighbors above which an (x,y) point
% is excluded

% xd and yd are the density-filtered points

% rho is the point-wise density (number of neighbors)


Mdl = createns([x y]);                      % Create nearest neighbor searcher object for kD-tree search (slightly faster than exhaustive pairwise comparisons)

neighbors = rangesearch(Mdl,[x y],n_dist);  % Find all neighbors w/in distance n_dist of each (x,y) point; note that each point is considered its own neighbor

rho = nan(length(x),1);                     % Stores number of neighbors for each point

for i = 1:length(x)
    if ~isnan(x(i))
        rho(i) = length(neighbors{i})-1;    % Subtract one because neighbors of point x,y include x,y
    end
end

% Filter out points

xd = x;
yd = y;

xd(rho >= n_thres) = NaN;
yd(rho >= n_thres) = NaN;

end
