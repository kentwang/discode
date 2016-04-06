function result = binologlik(Y, p)
	x = sum(Y);
	n = length(Y);
	result = gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1) + x * log(p) + (n - x) * log(1 - p);
end