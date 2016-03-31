function [X, Y] = bibetarnd(n, a, b, c)
% Generate n random sample from bivariate beta distribution
% Ingram Olkin and Ruixue Liu (2003)

% TODO change to randg of which the seed can be set.

    U = gamrnd(a, 1, n, 1);
    V = gamrnd(b, 1, n, 1);
    W = gamrnd(c, 1, n, 1);

    X = U * 1.0 ./ (U + W);
    Y = V * 1.0 ./ (V + W);
end