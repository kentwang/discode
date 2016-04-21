function result = sampleNuep2(Z, U, V, K_plus, L_plus, sigma_w)
% sample nuep from its posterior distribution with W integrated out
	[I, J] = size(Z);
	tz = Z(:);
	F = kron(V(:, 1:L_plus), U(:, 1:K_plus));

	mu = tz'*(eye(I*J) - F*inv(F'*F + eye(K_plus*L_plus)/sigma_w^2)*F')*ones(I*J, 1)/(I*J + 1);
	sigma = 1/(I*J + 1);
	result = normrnd(mu, sigma);
end