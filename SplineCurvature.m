function [zz,L0] = SplineCurvature(x0,y0)
% Measures linearity/tortuosity on an interpolant with fixed window length
% Modified code from James Hengenius
% 2017-08-29

%k = ~isnan(x0) & ~isnan(y0);

% Interp points equidistant along arc length
pt_density = 11.2; %per unit arcl length (assume cm for most trail data)
arcl = nansum( sqrt( diff(x0).^2 + diff(y0).^2 ) );
zz = interparc(roundn(pt_density*arcl,0),x0,y0,'linear');

L0 = nan(length(zz),1);

% Window size
res = 1/250; % Resolution scaling term; eg, window length = 2% total length
win = roundn(res*length(zz),0);

if mod(win,2)==0
    win = win+1; % Make window an odd length so it is centered on each point
end

% Loop over all positions
for i = 1+(win-1)/2 : length(L0)-(win-1)/2
%     if ~mod(i,100000);
%         disp(i./(length(L0)-(win-1)/2));
%     end
    x_w = zz(i-(win-1)/2 : i+(win-1)/2, 1) ;
    y_w = zz(i-(win-1)/2 : i+(win-1)/2, 2) ;
    
    ab = sqrt( (x_w(end) - x_w(1))^2 + (y_w(end) - y_w(1))^2 ); % Euclidean distance of window
    cd = nansum( sqrt( diff(x_w).^2 + diff(y_w).^2 ) ); % arc length of window
    L0(i) = ab/cd; %linearity of window
end



