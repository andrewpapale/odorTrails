function [ conc ] = odorant_C(m_x, m_y, s_x, s_y, r_xy , x, y )
%ODORANT_C generates continuous concentration values from a bivariate norm
%distribution with parameters:
%   m_x     - mean of X
%   m_y     - mean of Y
%   s_x     - stdev of X
%   s_y     - stdev of Y
%   r_xy    - correlation of X and Y
%   x       - x where function is evaluated
%   y       - y where function is evalutated
%   
%   

c_max = ( 2*pi()*s_x*s_y*sqrt( 1-r_xy ) ).^-1 * exp( -(2-2*r_xy.^2).^-1 * ( (0-m_x)^2/s_x^2 + (0-m_y)^2/s_y^2 - 2*r_xy*(0-m_x)*(0-m_y)/(s_x*s_y) )  );

conc  = ( 2*pi()*s_x*s_y*sqrt( 1-r_xy ) ).^-1 * exp( -(2-2*r_xy.^2).^-1 * ( (x-m_x)^2/s_x^2 + (y-m_y)^2/s_y^2 - 2*r_xy*(x-m_x)*(y-m_y)/(s_x*s_y) )  );
conc  = conc/c_max;
end

