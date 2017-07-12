function result = computeline(a,m,xr)

%m = (a(2) - b(2)) / (a(1) - b(1));
n = a(2) - a(1) * m;

x = xr(1):xr(2);
y = m * x + n;
y = round(y);

for i = 1 : length(x)
    result{i} = [x(i), y(i)];
end

end

