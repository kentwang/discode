function result = W_expand(new_k, dim, U, V, W, nuep, sigma_w)
% expand rows or columns of W according to the generation of
% new features of U or V. Individual element of the new vector of
% W is sampled.

% INPUT
%	dim - 2 for row and 1 for column
%	new_k - number of new rows or columns to append

% TEST
% W = W_true; U = U_true; V = V_true; % with new feature in sampleUV.m 

	[K, L] = size(W);
	% append rows
	if dim == 2
		W((K + 1):(K + new_k), :) = 0;
		% Gibbs  and MH sampling of new elements
		for d = 1:new_k
			for l = 1:L
				l_curr = log(normpdf(W(K + d, l), 0, sigma_w)) + Z_UVWlogpdf(U, V, W, nuep);
				w_prop = normrnd(W(K + d, l), sigma_w); % simply use normal proposal
				W(K + d, l) = w_prop;
				l_prop = log(normpdf(W(K + d, l), 0, sigma_w)) + Z_UVWlogpdf(U, V, W, nuep);
				acc_w = exp(min(0, l_prop - l_curr));

				% printf('iter d = %d, l = %d; l_curr = %f, l_prop = %f; w_prop = %f; accept ratio = %f \n', d, l, l_curr, l_prop, w_prop, acc_w);

				if rand < acc_w
					W(K + d, l) = w_prop;
				end
			end
		end
	elseif dim == 1
		W(:, (L + 1):(L + new_k)) = 0;
		for d = 1:new_k
			for k = 1:K
				l_curr = normpdf(W(k, L + d), 0, sigma_w) * exp(Z_UVWlogpdf(U, V, W, nuep));
				w_prop = normrnd(W(k, L + d), sigma_w);
				W(k, L + 1) = w_prop;
				l_prop = normpdf(W(k, L + d), 0, sigma_w) * exp(Z_UVWlogpdf(U, V, W, nuep));
				acc_w = exp(min(0, l_prop - l_curr));

				if rand < acc_w
					W(k, L + d) = w_prop;
				end
			end
		end
	end
	result = W;
end