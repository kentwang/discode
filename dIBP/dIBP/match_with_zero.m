function mu = match_with_zero(mu_u, mu_v)
% match the probability pairs with not really zero (singularity)
%
% TODO - fix the zero division problem of updateAB.m

	n_u = length(mu_u);
	n_v = length(mu_v);

	if n_u > n_v
		for i = (n_v + 1):(n_u)
			mu_v = [mu_v; mu_v(i - 1) / 3];
		end
		% mu_v = [mu_v; .00001 * ones(n_u - n_v, 1)];
	elseif n_u < n_v
		for i = (n_u + 1):(n_v)
			mu_u = [mu_u; mu_u(i - 1) / 3];
		end
		% mu_u = [mu_u; .00001 * ones(n_v - n_u, 1)];
	end
	mu = [mu_u, mu_v];
end