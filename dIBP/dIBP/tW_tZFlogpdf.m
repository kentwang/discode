function result = tW_tZFlogpdf(W, Z, U, V, nuep, sigma_w)
% log pdf of posterior of stacked W. May be able to get around using cancelation.

% NOTE
%	- |SIGMA| will explode is sigma_w is too large!!!

	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tz = Z(:);

	SIGMA = inv(F' * F + eye(K * L, K * L) / sigma_w^2);
	
	MU = SIGMA * F' * (tz - nuep * ones(I * J, 1));

	result = log(mvnpdf(W(:)', MU', SIGMA));
end