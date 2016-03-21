function [a, b] = sampleAB(I, J, K, L)
% update alpha and beta
	a = mygamrnd(1 + K, 1/(1 + harmonic(I)), 1);
	b = mygamrnd(1 + L, 1/(1 + harmonic(J)), 1);
end