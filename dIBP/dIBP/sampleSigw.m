function [sigma_w, a_sigw, b_sigw] = sampleSigw(a_sigw, b_sigw, K, L, W)
% sample sigma_w from the posterior distribution

% NOTE
%	- gamrnd use a, b for shape and scale. We keep using scale (1/beta)
%	- but invgamma-normal conjugacy use shape and rate

	a_sigw = a_sigw + K * L / 2;
	b_sigw = 2 * b_sigmw / (2 + b_sigmw * prod(W(:)));
	sigma_w = gamrnd(a_sigw, b_sigmw);
end