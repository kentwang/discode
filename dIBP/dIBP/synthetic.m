function [U_orig, V_orig, Z_orig, X] = synthetic(W, I, J, K, L, nuep)
% Simulate ground truth data
% Output:
%   - U, V,
%   - Z, X
    
    % U = featrnd(I, K, 0.7);
    % V = featrnd(J, L, 0.8);

    U_orig = zeros(I, K);
    V_orig = zeros(J, L);

	for i = 1:I
		U_orig(i, :) = (rand(1, K) > .5);
		while(sum(U_orig(i, :)) == 0)
			U_orig(i, :) = (rand(1, K) > .5);
		end
	end

    for j = 1:J
    	V_orig(j, :) = (rand(1, L) > .5);
    	while(sum(V_orig(j, :)) == 0)
    		V_orig(j, :) = (rand(1, L) > .5);
    	end
    end

    Eps = normrnd(nuep, 1, I, J);
    Z_orig = U_orig * W * V_orig' + Eps;
    X = (Z_orig > 0);
end