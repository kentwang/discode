function output = UV_mupdf(U, V, mu_u, mu_v)
% ProbUV_mu - compute the probability Pr(U, V | \mu) in Eq. (2.20)
% Note: mu_u and mu_v are of the real active #features
    [I, ~] = size(U);
    [J, ~] = size(V);
    
    output = prod(mu_u.^(sum(U, 1)')) * ...
        prod((1 - mu_u).^(I - sum(U, 1)')) * ...
        prod(mu_v.^(sum(V, 1)')) * ...
        prod((1 - mu_v).^(J - sum(V, 1)'));
end

