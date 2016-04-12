function output = Z_UVlogpdf(Z, U, V, sigma_w, nuep, I, J, K_plus, L_plus)

% TODO
%	- Might be improved by efficient updating of Kronecker product

	tz = Z(:);
	F = kron(V, U);
	SIGMA = inv(F'*F + eye(K_plus*L_plus, K_plus*L_plus));

	output = -.5*I*J*log(2*pi) - K_plus*L_plus*log(sigma_w) - .5*logdet(SIGMA) - .5*(tz-nuep*ones(I*J, 1))'*(eye(I*J)-F*SIGMA*F')*(tz-nuep*ones(I*J, 1));
end