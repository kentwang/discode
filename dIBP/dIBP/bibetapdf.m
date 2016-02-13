% Test
% a = 2; b = 2; c = 2;
% [x, y] = meshgrid([0.001:0.01:0.999]);
% Z = bibetapdf(x, y, a, b, c);
% surf(x, y, Z, gradient(Z));
% colorbar;

function output = bibetapdf(x, y, a, b, c)
% bi-variate beta pdf in Olkin and Liu [2003]
    output = x.^(a - 1) .* y.^(b - 1) .* (1 - x).^(b + c - 1) .* (1 - y).^ ...
        (a + c - 1) / gamma(a) / gamma(b) / gamma(c) * gamma(a + b + c) ./ ...
        (1 - x .* y).^(a + b + c); 
end

