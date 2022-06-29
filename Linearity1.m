function L1 = Linearity1(t,y,window)
nW = length(t);
L1 = nan(nW,1);
for iW=1:(nW-window)
    x_w = t(iW:iW+window);
    y_w = y(iW:iW+window);
    ab = sqrt((x_w(end) - x_w(1))^2 + (y_w(end) - y_w(1))^2 );
    %cd = nansum(sqrt( diff(x_w).^2 + diff(y_w).^2 ) );
    cd = nansum(sqrt(x_w.^2 + y_w.^2));
    L1(iW) = ab/cd;
end

end

