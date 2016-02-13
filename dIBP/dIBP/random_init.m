rng(20160209);

I = 10;
J = 20;
K = 4;
L = 6;
M = max(K, L);

a = 2; b = 5; c = 1;
[lambda_u, lambda_v] = bibetarnd(M, a, b, 1);
mu_u = cumprod(lambda_u);
mu_v = cumprod(lambda_v);
mu = [mu_u, mu_v];

mu_u_test = mu_u + normrnd(0, 0.001, M, 1);
mu_v_test = mu_v + normrnd(0, 0.001, M, 1);
mu_test = [mu_u_test, mu_v_test];

U = matbinornd(mu_u, I, K);
V = matbinornd(mu_v, J, L);

UV_mupdf(U, V, mu_u(1:K), mu_v(1:L)) / ...
    UV_mupdf(U, V, mu_u_test(1:K), mu_v_test(1:L))

%-- Test computation of acceptance ratio
UV_mupdf(U, V, mu_u_test(1:K), mu_v_test(1:L)) / ...
    UV_mupdf(U, V, mu_u(1:K), mu_v(1:L)) * ...
    mucondpdf(mur, murp, murm, a, b, c)



