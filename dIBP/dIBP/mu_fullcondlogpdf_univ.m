function result = mu_fullcondlogpdf_univ(mu_1, mu_2, mu_3, a)
% the conditional log pdf of mu_1 | mu_2, mu_3 (univariate version)
	result = log(betapdf(mu_2/mu_1, a, 1)) - log(mu_1) + log(betapdf(mu_3/mu_2, a, 1)) - log(mu_2);
end
