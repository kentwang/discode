function result = tW_tZFpdf(W, Z, U, V, nuep, sigma_w)
% pdf of posterior of stacked W

	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tz = Z(:);


	MU = inv(F' * F + eye(K * L, K * L) / sigma_w^2) * F' * (tz - nuep * ones(I * J, 1));
	SIGMA = inv(F' * F + eye(K * L, K * L) / sigma_w^2);

	result = mvnpdf(W(:)', MU', SIGMA);
end