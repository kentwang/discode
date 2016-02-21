function output = truncbetarnd(a, b, lower, upper)
% Simulate data from truncated beta distribution using CDF inverse
% Madarajah & Kotz 2006

% Todo
%    - test the random number generator

    p = unifrnd(0, 1);
    output = betaincinv(betacdf(lower, a, b) + p * (betacdf(upper, a, b) - ...
        betacdf(lower, a, b)), a, b);
    
end