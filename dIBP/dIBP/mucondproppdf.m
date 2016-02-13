function output = mucondproppdf(mur, murp, murm, a, b, c, K)
% Proposal pdf with truncated beta's in Eq (2.21).
% Note: check the equations
    
    output = truncbetapdf(mur(1), a/K, c, murp(1), murm(1)) * ...
        truncbetapdf(mur(2), b/K, c, murp(2), murm(2));
end
