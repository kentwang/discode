function [a, b] = updateAB(a, b, mu_u, mu_v)
% update alpha and beta using MH sampling with uniform perturbation

% NOTE
% 	- Affected by matched pairs of bibeta and unmatched univariate beta
	
	K = length(mu_u);
	L = length(mu_v);
	M_low = min(K, L);
    M_high = max(K, L);
    mu = [mu_u(1:M_low), mu_v(1:M_low)];

	% printf('[a, b]: [%f, %f]\n', a, b);
	% propose a and b using uniform perturbation
	if rand < .5
		a_prop = a + rand/20;
	else
		a_prop = a - rand/20;
	end

	if rand < .5
		b_prop = b + rand/20;
	else
		b_prop = b - rand/20;
	end

	if K < M_low
		logl_diff = log(mujoinpdf(mu, a_prop, b_prop)) + sum(log(betapdf(mu_v(M_low:M_high), b_prop, 1))) + log(gampdf(a_prop, 1, 1)) + log(gampdf(b_prop, 1, 1)) - log(mujoinpdf(mu, a, b)) - sum(log(betapdf(mu_v(M_low:M_high), b, 1))) - log(gampdf(a, 1, 1)) - log(gampdf(b, 1, 1));
	elseif L < M_low
		logl_diff = log(mujoinpdf(mu, a_prop, b_prop)) + sum(log(betapdf(mu_u(M_low:M_high), a_prop, 1))) + log(gampdf(a_prop, 1, 1)) + log(gampdf(b_prop, 1, 1)) - log(mujoinpdf(mu, a, b)) - sum(log(betapdf(mu_u(M_low:M_high), a, 1))) - log(gampdf(a, 1, 1)) - log(gampdf(b, 1, 1));
	else
		logl_diff = log(mujoinpdf(mu, a_prop, b_prop)) + log(gampdf(a_prop, 1, 1)) + log(gampdf(b_prop, 1, 1)) - log(mujoinpdf(mu, a, b)) - log(gampdf(a, 1, 1)) - log(gampdf(b, 1, 1));
	end

	ab_acc = exp(min(0, logl_diff));

	if rand < ab_acc
		a = a_prop;
		b = b_prop;
	end
end