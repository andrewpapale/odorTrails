function dx = logError(x,dim)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

assert(dim==1 | dim==2,'error dimension must be 1 or 2');

if dim==1
    dim1=2;
elseif dim==2
    dim1=1;
end

dx = log10(1+(nanstd(x,[],dim)./(nanmean(x,dim).*sqrt(size(x,dim1)))));


end

