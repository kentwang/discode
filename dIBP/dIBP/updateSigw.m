function result = updateSigw(sigma_w, U, V, W, Z, nuep)
% update sigma_w using MH sampling

% TODO
% 	- Cancel out common terms when propose and sample.
	
	% propose sigma_w using uniform perturbation
	if rand < .5
		sigma_w_prop = sigma_w + rand/20;
	else
		sigma_w_prop = sigma_w - rand/20;
	end

	% sigma_w_acc = min(1, tW_tZFpdf(W, Z, U, V, nuep, sigma_w_prop) / tW_tZFpdf(W, Z, U, V, nuep, sigma_w));
	sigma_w_acc = exp(min(0, tW_tZFlogpdf(W, Z, U, V, nuep, sigma_w_prop) - tW_tZFlogpdf(W, Z, U, V, nuep, sigma_w)));

	if rand < sigma_w_acc
		result = sigma_w_prop;
	else
		result = sigma_w;
	end
end