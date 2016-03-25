function result = mujoinpdf(mu, a, b)
% Joint pdf of a pair of vector mu_u and mu_v. Return a real-valued density
% This function is used for sampling from posterior of a and b

% Needs to be checked and tested

    n = size(mu, 1);
    if n == 1
        result = bibetapdf(mu(1), mu(2), a, b, 1)
    else
        result = zeros(n, 1);
        result(1) = bibetapdf(mu(1, 1), mu(1, 2), a, b, 1);
        for i=2:n
            difmu = mu(i, :) ./ mu(i-1, :);
            result(i) = bibetapdf(difmu(1), difmu(2), a, b, 1) / prod(mu(i-1, :));
            % printf('value of p_i: %f \n', result(i));
            % printf('values of diff diff mu_i: [%f, %f] \n', difmu(1), difmu(2));
        end
        result = prod(result);
    end
end