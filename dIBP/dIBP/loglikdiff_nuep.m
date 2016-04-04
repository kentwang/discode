function result = loglikdiff_nuep(W, Z, U, V, sigma_w, nuep, nuep_prop)
% calculate the difference in log likelihood by sigma_w
	
% TODO
% 	- DO I REALLY NEED THIS?
%	- sigma_w has a close-form posterior distribution.

	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tz = Z(:);
	bOne = ones(I * J, 1);

	SIGMA = inv(F' * F + eye(K * L, K * L) / sigma_w^2);

	result = (nuep_prop - nuep) * (((tz' * F) * SIGMA) * (F' * bOne) - W(:)' * (F' * bOne)) - ...
	.5 * (nuep_prop^2 - nuep^2) * ((bOne' * F) * SIGMA) * (F' * bOne);
end	