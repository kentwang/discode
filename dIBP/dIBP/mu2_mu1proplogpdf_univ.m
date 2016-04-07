function result = mu2_mu1proplogpdf_univ(mu_1, mu_2, a)
% proposed conditional pdf of the LAST pair of mu_u or mu_v

% NOTE
%	- a and b should be a/K and b/L in real sampling. See Xuan[2015]

	result = truncbetalogpdf(mu_2, a, 1, 0, mu_1);
end