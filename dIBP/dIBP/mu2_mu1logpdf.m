function result = mu2_mu1logpdf(mu_1, mu_2, a, b)
% conditional pdf of the LAST pair of mu
	
	result = bibetalogpdf(mu_2(1)/mu_1(1), mu_2(2)/mu_1(2), a, b, 1) - sum(log(mu_1));
end