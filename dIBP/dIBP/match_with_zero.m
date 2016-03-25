function mu = match_with_zero(mu_u, mu_v)
% match the probability pairs with not really zero (singularity)

	n_u = length(mu_u);
	n_v = length(mu_v);

	if n_u >= n_v
		mu_v = [mu_v; .00001 * ones(n_u - n_v, 1)];
	else
		mu_u = [mu_u; .00001 * ones(n_v - n_u, 1)];
	end
	mu = [mu_u, mu_v];
end