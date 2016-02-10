function U = matbinornd(mu_u, I, K)
% Generate binary feature matrix
% Note: mu_u is already cumulative prod and full
    U = zeros(I, K);
    for k = 1:K
        U(:, k) = binornd(1, mu_u(k), I, 1);
    end
end