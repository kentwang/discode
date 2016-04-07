function output = truncbetalogpdf(x, a, b, lower, upper)
% Truncated beta pdf
    
    C = betacdf(upper, a, b) - betacdf(lower, a, b);
    output = log(betapdf(x, a, b) / C);
end