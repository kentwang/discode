function output = truncbetapdf(x, a, b, lower, upper)
% Truncated beta pdf
    
    C = betacdf(upper, a, b) - betacdf(lower, a, b);
    output = betapdf(x, a, b) / C;
end