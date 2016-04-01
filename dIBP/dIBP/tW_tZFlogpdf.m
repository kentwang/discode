function result = tW_tZFlogpdf(W, Z, U, V, nuep, sigma_w)
% pdf of posterior of stacked W

% NOTE
%	- |SIGMA| will explode is sigma_w is too large!!!

	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tz = Z(:);


	MU = inv(F' * F + eye(K * L, K * L) / sigma_w^2) * F' * (tz - nuep * ones(I * J, 1));
	SIGMA = inv(F' * F + eye(K * L, K * L) / sigma_w^2);

	result = log(mvnpdf(W(:)', MU', SIGMA));
end