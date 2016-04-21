function result = sampleSigw2(Z, U, V, K_plus, L_plus, sigma_w, nuep)
% TODO
%	- use log version of mvnpdf

	[I, J] = size(Z);
	tz = Z(:);
	F = kron(V(:, 1:L_plus), U(:, 1:K_plus));

	MU = nuep * ones(I*J, 1);
	SIGMA_curr = inv(eye(I*J) - F*inv(F'*F + eye(K_plus*L_plus)/sigma_w^2)*F');
    loglik_curr = logmvnpdf(tz', MU', SIGMA_curr);


	% propose sigma_w using uniform perturbation
	if rand < .5
		sigma_w_prop = sigma_w + rand/20;
	else
		sigma_w_prop = sigma_w - rand/20;
	end

	SIGMA_prop = inv(eye(I*J) - F*inv(F'*F + eye(K_plus*L_plus)/sigma_w_prop^2)*F');
	loglik_prop = logmvnpdf(tz', MU', SIGMA_prop);

	sigma_w_acc = exp(min(0, loglik_prop - loglik_curr));

	if rand < sigma_w_acc
		result = sigma_w_prop;
	else
		result = sigma_w;
	end
end