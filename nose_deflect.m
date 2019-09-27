function [ th_head, th_nose, th_bn ] = nose_deflect( x0, y0, nx0, ny0)
%nose_deflect takes x,y trajectories for nose nx0,ny0 and body x0,y0. For
% frame indices t=2:end, it computes an estimate of forward heading th_head, the
%relative nose direction (th_nose) and the angular difference between the
%two th_bn.

%th_bn is the measure I have been using.


%% What's being calculated
% The forward heading of the mouse center of mass (x0,y0) is given as:

%       th_head(t) = atan2( y0(t) - y0(t-1) ,  x0(t) - x0(t-1)   ),

% the movement vector from the previous to current frame.
% Note: This is noisy. Using the mean coordinates of the last n framws can
% help stabilize it:

%       th_head(t) = atan2( y0(t) - nanmean(y0(t-n:t-1)) ,  x0(t) - nanmean(x0(t-n:t-1))   ),



% The nose orientation of the mouse is given as:
%       th_nose(t) = atan2( ny0(t) - y0(t) ,  nx0(t) - x0(t)   )



% The nose deflection relative to forward heading is given as:
%       th_bn(t) = angdiff(   th_head(t) , th_nose(t)   )


%% Check for unequal vector lengths
if length(x0) ~= length(y0) | length(x0) ~= length(nx0) | length(x0) ~= length(ny0)
    error('Body and nose coordinate vectors must be the same length')
end

%% Compute values to return
th_head = nan(length(x0),1);
th_nose = nan(length(x0),1);
th_bn   = nan(length(x0),1);

%For smoother estimate of th_head, increase dt_gap
dt_gap = 10; % number of previous frames to average over when estimating forward heading

for t = dt_gap+1:length(x0)
%     th_head(t) =    atan2( y0(t) - y0(t-1) ,  x0(t) - x0(t-1)   );
    th_head(t) =    atan2d(  y0(t)-nanmean( y0(t-dt_gap:t-1) )     , x0(t)-nanmean(x0(t-dt_gap:t-1))   );
    th_nose(t) =    atan2d( ny0(t) - y0(t) ,  nx0(t) - x0(t)   );
    th_bn(t)   =    th_nose(t)-th_head(t);
end

