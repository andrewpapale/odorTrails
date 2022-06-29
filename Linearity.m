function [L,AB,CD] = Linearity(x,y)
% 2017-08-15 AndyP
% [L,seglen] = Linearity(x,y,window);
% linear distance / summed distance along arc length

nW = length(x);
AB = nan(size(x));
% compute velocity using adaptive windowing method, smooth at 0.25s
[dx,winusedx] = foaw_diff(x,1/50,50,0.5,0.25);
[dy,winusedy] = foaw_diff(y,1/50,50,0.5,0.25);
CD = sqrt(1+(dy./dx).^2); % compute arc length at each point
for iW=1:nW
    x1 = x(iW);
    y1 = y(iW);
    if ~isnan(winusedx(iW)) & ~isnan(winusedy(iW)) %#ok<AND2>
        x2 = x(min(iW+winusedx(iW),length(x)));
        y2 = y(min(iW+winusedy(iW),length(y)));
        AB(iW) = sqrt((x2-x1).^2+(y2-y1).^2); % compute "linear" distance
    end
end
L = AB./CD;
end

