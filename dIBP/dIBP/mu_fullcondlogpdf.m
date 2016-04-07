function result = mu_fullcondlogpdf(mu_1, mu_2, mu_3, a, b)
% the conditional log pdf of mu_1 | mu_2, mu_3
	result = bibetalogpdf(mu_2(1)/mu_1(1), mu_2(2)/mu_1(2), a, b, 1) - sum(log(mu_1)) + bibetalogpdf(mu_3(1)/mu_2(1), mu_3(2)/mu_2(2), a, b, 1) - sum(log(mu_2));
end
