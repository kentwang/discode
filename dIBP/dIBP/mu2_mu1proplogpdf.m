function result = mu2_mu1proplogpdf(mu_1, mu_2, a, b)
% proposed conditional pdf of the LAST pair of mu

% NOTE
%	- a and b should be a/K and b/L in real sampling. See Xuan[2015]

	result = log(truncbetapdf(mu_2(1), a, b, 0, mu_1(1))) + log(truncbetapdf(mu_2(2), a, b, 0, mu_1(2)));
end