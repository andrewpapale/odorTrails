function [c] = odorFun_grid(x,y,cgrid)


side_lx = 114.3/2;
side_ly = 91.400/2;


x_vals = -side_lx:0.1:side_lx;
y_vals = -side_ly:0.1:side_ly;

[~,x_idx] = min(abs(x_vals - x));
[~,y_idx] = min(abs(y_vals - y));


    
    c = cgrid(x_idx,y_idx);