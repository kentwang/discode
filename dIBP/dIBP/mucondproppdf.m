function output = mucondproppdf(mur, murp, murm, a, b, c, K_plus, L_plus)
% Proposal pdf with truncated beta's in Eq (2.21).
% Note: check the equations

% Todo:
%   - verify the consistency of -1

    if ~any(mur == -1)
        % both muru and murv exist and two truncated beta are used
        output = truncbetapdf(mur(1), a/K_plus, c, murp(1), murm(1)) * ...
            truncbetapdf(mur(2), b/L_plus, c, murp(2), murm(2));
    elseif find(mur == -1) == 1
        % muru is null and truncbeta(b/L, c) is used
        output = truncbetapdf(mur(2), b/L_plus, c, murp(2), murm(2));
    elseif find(mur == -1) == 2
        % murv is null and truncbeta(a/K, c) is used
         output = truncbetapdf(mur(1), a/K_plus, c, murp(1), murm(1));
    end
end
