function output = mujoinpdf(mu, a, b, c)
% Joint pdf of mu_u and mu_v. Eq. (19) second equality sign
% Output: return a density vector of length of mu

% Note: THIS FUNCTION IS NOT NECESSARY!

    n = size(mu, 1);
    if n == 1
        output = bibetapdf(mu(1), mu(2), a, b, c);
    else
        output = zeros(n, 1);
        output(1) = bibetapdf(mu(1, 1), mu(1, 2), a, b, c);
        for i=2:n
            difmu = mu(i, :) ./ mu(i-1, :);
            output(i) = bibetapdf(difmu(1), difmu(2), a, b, c) / prod(mu(i-1, :));
        end
    end
end