function result = mu1_mu2proplogpdf(mu_1, mu_2, a, b)
% proposed conditional pdf of the first pair of mu

% NOTE
%	- a and b should be a/K and b/L in real sampling. See Xuan[2015]

	result = truncbetalogpdf(mu_1(1), a, b, mu_2(1), 1) + truncbetalogpdf(mu_1(2), a, b, mu_2(2), 1);
end