function result = sampleZ(X, U, V, W, nuep)
% sample Z using its posterior distirbution.
% If you wish to simulate a random variable
% 'Z' from the non-standard Gaussian N(m,s^2)
% conditional on l<Z<u, then first simulate
% X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;

	
	[I, J] = size(X);
	result = zeros(I, J);

	for i = 1:I
		for j = 1:J
			z_mu = U(i, :) * W * V(j, :)' + nuep;
			z_sig = 1;
			if X(i, j) == 1
				lower = 0;
				upper = Inf;
				result(i, j) = z_mu + z_sig * trandn((lower - z_mu) / z_sig, (upper - z_mu) / z_sig);
			else
				lower = -Inf;
				upper = 0;
				result(i, j) = z_mu + z_sig * trandn((lower - z_mu) / z_sig, (upper - z_mu) / z_sig);
			end	
		end
	end
end