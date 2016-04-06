function result = mu1_mu2logpdf(mu_1, mu_2, a, b)
% conditional pdf of the first pair of mu
	result = bibetalogpdf(mu_2(1)/mu_1(1), mu_2(2)/mu_1(2), a, b, 1) + sum(log(mu_2)) - 2*sum(log(mu_1));
end