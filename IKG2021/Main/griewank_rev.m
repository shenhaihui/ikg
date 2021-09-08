function y = griewank_rev(i,x)

[d, m]=size(x); % d - dimension, m - number of points
y_temp = 1.5^(d-1) * prod(cos(x ./ repmat(sqrt(i*(1:d)'),1,m)),1);
y = sum(x.^2,1)/4000 - y_temp;

end