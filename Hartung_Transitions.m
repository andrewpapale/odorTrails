% 2019-10-10 AndyP
% get transitions for Hartung data + maze

k0 = confb > 0.3 & ~(Top & Right) & ~(Top & Left) & ~(Bottom & Right) & ~(Bottom & Left) & ~(~Bottom & ~Top & ~Left & ~Right & ~center) & ~isinf(log10dphi) & sqrt((nx-x).^2+(ny-y).^2) < 80;

center = (x-nanmedian(x)).^2+(y-nanmedian(y)).^2 < 700;
Top = y > 320 & ~center;
Bottom = y < 280 & ~center;
Left = x < 360 & ~center;
Right = x > 400 & ~center;



for iS=1:22
    k = sess==iS & k0;
    
    T0 = Top(k);
    B0 = Bottom(k);
    L0 = Left(k);
    R0 = Right(k);
    x0 = x(k);
    y0 = y(k);
    nx0 = nx(k);
    ny0 = ny(k);
    c0 = center(k);
    
    
    
    
    

end
        
    
    
