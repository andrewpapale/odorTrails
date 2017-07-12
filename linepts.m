function result = linepts(a, b)

m = (a(2) - b(2)) / (a(1) - b(1));
n = b(2) - b(1) * m;

x = min(a(1), b(1)) : max(a(1), b(1));
y = m * x + n;
y = round(y);

for i = 1 : length(x)
    result{i} = [x(i), y(i)];
end

end

