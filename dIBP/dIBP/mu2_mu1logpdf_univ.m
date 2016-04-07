function result = mu2_mu1logpdf_univ(mu_1, mu_2, a)
% conditional pdf of the LAST pair of mu_u or u_v
	
	result = log(betapdf(mu_2(1)/mu_1(1), a, 1)) - log(mu_1);
end