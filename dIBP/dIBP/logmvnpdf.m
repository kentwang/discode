function result = logmvnpdf(x, mu, sigma)
	p = length(x);
	r = chol(sigma);
	result = -p*log(2*pi)/2 - sumsq((x-mu)/r, 2)/2 - sum(log(diag(r)));
end