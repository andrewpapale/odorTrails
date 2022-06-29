 
ix0 = 0.001:0.001:1;
is0 = 0.001:0.001:0.5;
for ix=1:length(ix0)
    for is=1:length(is0) 
        G(ix,is) = exp(-ix0(ix)./(2*(is0(is)).^2)); 
        P(ix,is) = ix0(ix).^(-is0(is));
        L(ix,is) = -is0(is).*ix0(ix);
    end 
end