%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing sampling algorithm using random feature matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20160209);

I = 10;
J = 20;
K = 4;
L = 6;
M = max(K, L);

% sample the feature generating probabilities with augmentation
a = 2; b = 5; c = 1;
[lambda_u, lambda_v] = bibetarnd(M, a, b, 1);
mu_u_aug = cumprod(lambda_u);
mu_v_aug = cumprod(lambda_v);
mu = [mu_u_aug, mu_v_aug];
mu_u = mu_u_aug(1:K);
mu_v = mu_v_aug(1:L);

mu_u_test_aug = mu_u_aug + normrnd(0, 0.001, M, 1);
mu_v_test_aug = mu_v_aug + normrnd(0, 0.001, M, 1);
mu_test = [mu_u_test_aug, mu_v_test_aug];
mu_u_test = mu_u_test_aug(1:K);
mu_v_test = mu_v_test_aug(1:L);


U = matbinornd(mu_u, I, K);
V = matbinornd(mu_v, J, L);

%-- Test computation of acceptance ratio for paired case
r = 2;
mur = mu(r, :);
murp = mu(r + 1, :);
murm = mu(r - 1, :);
mur_test = mu_test(r, :);
murp_test = mu_test(r + 1, :);
murm_test = mu_test(r - 1, :);

UV_mupdf(U(:, r), V(:, r), mu_u_test(r), mu_v_test(r)) / ...
    UV_mupdf(U(:, r), V(:, r), mu_u(r), mu_v(r)) * ...
    mucondpdf(mur_test, murp_test, murm_test, a, b, c) / ...
    mucondpdf(mur, murp, murm, a, b, c) * ...
    mucondproppdf(mur, murp, murm, a, b, c, K, L) / ...
    mucondproppdf(mur_test, murp_test, murm_test, a, b, c, K, L)

%-- Test computation of acceptance ratio for unpaired case
%-- U probability is nullified
mur(1) = -1; 
murp(1) = -1; 
murm(1) = -1;
mur_test(1) = -1; 
murp_test(1) = -1; 
murm_test(1) = -1;
    
UV_mupdf(U(:, r), V(:, r), mu_u_test(r), mu_v_test(r)) / ...
    UV_mupdf(U(:, r), V(:, r), mu_u(r), mu_v(r)) * ...
    mucondpdf(mur_test, murp_test, murm_test, a, b, c) / ...
    mucondpdf(mur, murp, murm, a, b, c) * ...
    mucondproppdf(mur, murp, murm, a, b, c, K, L) / ...
    mucondproppdf(mur_test, murp_test, murm_test, a, b, c, K, L)


