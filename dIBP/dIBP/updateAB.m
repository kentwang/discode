function [a, b] = updateAB(a, b, mu)
% update alpha and beta using MH sampling with uniform perturbation
	
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

	ab_acc = mujoinpdf(mu, a_prop, b_prop) * gampdf(a_prop, 1, 1) * gampdf(b_prop, 1, 1) ...
	/ mujoinpdf(mu, a, b) / gampdf(a, 1, 1) / gampdf(b, 1, 1);

	if rand < ab_acc
		a = a_prop;
		b = b_prop;
	end
end