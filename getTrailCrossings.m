function tc = getTrailCrossings(x0,y0,xT,yT)
% 2017-07-19 AndyP 
% tc = getTrailCrossings(x0,y0,xT,yT)
% function to get trail crossing times

tol = 1;

nP = length(x0);
tc = zeros(nP,1);
for iP=1:nP
    for iT=1:length(xT)
        if sqrt((x0(iP)-yT(iT)).^2 + (y0(iP)-xT(iT)).^2)<tol
            tc(iP)=1;
        end
    end
end

end

