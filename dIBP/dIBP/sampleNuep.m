function result = samplenuep(U, V, W, Z)
% sample nuep directly from the posterior distribution

	[I, J] = size(Z);
	tw = W(:);
	tz = Z(:);
	F = kron(V, U);

	MU = (tz' * ones(I * J, 1) - ones(I * J, 1)' * F * tw) / (I * J + 1);
	SIGMA = (I * J + 1)^(-.5); 

	result = mvnrnd(MU', SIGMA);
end