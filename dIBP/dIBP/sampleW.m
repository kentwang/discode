function result = sampleW(Z, U, V, nuep, sigma_w)
% sample W from its stacked posterior distribution

% TEST
% U = U_true; V = V_true; Z = Z_true;

	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tZ = Z(:);


	SIGMA = inv(F' * F + eye(K * L, K * L) / sigma_w^2);
	MU = SIGMA * F' * (tZ - nuep * ones(I * J, 1));
	tW_new = mvnrnd(MU, SIGMA)';

	result = reshape(tW_new, K, L);
end
