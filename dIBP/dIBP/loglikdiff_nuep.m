function result = loglikdiff_nuep(W, Z, U, V, sigma_w, nuep, nuep_prop)
% calculate the difference in log likelihood by sigma_w
	
	[I, J] = size(Z);
	K = size(U, 2); L = size(V, 2);
	F = kron(V, U);
	tz = Z(:);
	bOne = ones(I * J, 1);
	bEye = eye(K * L, K * L);

	


end	