function result = sampleZ2(X, U, V, K_plus, L_plus, sigma_w, nuep)
% sample Z using its posterior distirbution with W integrated out

% If you wish to simulate a random variable
% 'Z' from the non-standard Gaussian N(m,s^2)
% conditional on l<Z<u, then first simulate
% X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;

% TODO
%	- sigma_w updating might be added

	tx = X(:);	
	[I, J] = size(X);
	result = zeros(I, J);

	F = kron(V(:, 1:L_plus), U(:, 1:K_plus));
	MU = nuep * ones(length(tx), 1);
	SIGMA = inv(eye(I*J) - F*(F'*F + eye(K_plus*L_plus)/sigma_w^2)^-1*F');
	diagSIGMA = diag(SIGMA);

	
	lower = zeros(length(tx), 1);
	lower(tx == 0) = -Inf;
	upper = zeros(length(tx), 1);
	upper(tx == 1) = Inf;

	result = MU + diagSIGMA.*trandn((lower-MU)./diagSIGMA, (upper-MU)./diagSIGMA);
end