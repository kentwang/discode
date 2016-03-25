function [U, V, W] = sampleUV(Z, U, V, W, nuep, a, b, sigma_w)
% Update U and V (simultaneously K and L)
% Augment W for new features of U and V

% Test
% U = U_true; V = V_true; W = W_true; Z = Z_true;
    
    % update U (and W if any new row features from U)
    [I, K] = size(U);
    for i = 1:I
        for k = 1:K
            if k > K
                break;
            end
            
            % shrink singular column
            % Note: if all zeros except U(i, k), then p(u_{ik}=1|) = 0!
            % Todo: - what if k = K  and singular?
            if U(i, k) > 0
                if U(:, k) - U(i, k) == 0
                    U(i, k) = 0;
                    U(:, k:(K-1)) = U(:, (k + 1):K);
                    K = K - 1;
                end
            end
                    
            % sample U
            U(i, k) = 1;
            p(1) = Z_UVWlogpdf(Z, U, V, W, nuep) + log(sum(U(:, k) - 1)) - ...
                log(I);
            U(i, k) = 0;
            p(2) = Z_UVWlogpdf(Z, U, V, W, nuep) + log(I - sum(U(:, k))) - ...
                log(I);
            p = exp(p - max(p));
            
            if rand < p(1)/(p(1)+p(2))
                U(i,k) = 1;
            else
                U(i,k) = 0;
            end
        end

        % sample new features of U
        n_new_lim = 4;
        trunc = zeros(1, n_new_lim + 1);
        W_hist = {};
        alpha_I = a / I; % alpha/I
        
        for k_i = 0:n_new_lim
            U(i, (K + 1):(K + k_i)) = 1; % test here U with new feature
            W_aug = W_expand(k_i, 2, Z, U, V, W, nuep, sigma_w); % append rows for new features U. TEST GOOD!
            trunc(k_i + 1) = k_i * log(alpha_I) - alpha_I - ...
                log(factorial(k_i)) + Z_UVWlogpdf(Z, U, V, W_aug, nuep);
            W_hist{k_i + 1} = W_aug;
        end
        
        U(i, (K + 1):(K + n_new_lim)) = 0;
        trunc = exp(trunc - max(trunc));
        trunc = trunc/sum(trunc);
        p = rand;
        t = 0;
        for k_i=0:n_new_lim
            t = t+trunc(k_i+1);
            if p < t
                new_dishes = k_i;
                break;
            end;
        end;
        U(i, (K + 1):(K + new_dishes)) = 1; % exchangeable
        U = U(:, 1:(K + new_dishes)); % shink the zero columns
        W = W_hist{new_dishes + 1};
        K = K + new_dishes;
    end

    % plot the updated U matrix
    % clf; imagesc(U);

    % update V (and W if any new column features from V)
    [J, L] = size(V);
    for j = 1:J
        for l = 1:L
            if l > L
                break;
            end

            if V(j, l) > 0
                if V(:, l) - V(j, l) == 0
                    V(j, l) = 0;
                    V(:, l:(L-1)) = V(:, (l + 1):L);
                    L = L - 1;
                end
            end

            % sample V
            V(j, l) = 1;
            p(1) = Z_UVWlogpdf(Z, U, V, W, nuep) + log(sum(V(:, l) - 1)) - ...
                log(J);
            V(j, l) = 0;
            p(2) = Z_UVWlogpdf(Z, U, V, W, nuep) + log(J - sum(V(:, l))) - ...
                log(J);
            p = exp(p - max(p));
            
            if rand < p(1)/(p(1)+p(2))
                V(j,l) = 1;
            else
                V(j,l) = 0;
            end
        end

        % samlpe new features of V
        n_new_lim = 4;
        trunc = zeros(1, n_new_lim + 1);
        W_hist = {};
        beta_J = b / J; % beta/J
        
        for l_j = 0:n_new_lim
            V(j, (L + 1):(L + l_j)) = 1; % test here V with new feature
            W_aug = W_expand(l_j, 1, Z, U, V, W, nuep, sigma_w);
            trunc(l_j + 1) = l_j * log(beta_J) - beta_J - ...
                log(factorial(l_j)) + Z_UVWlogpdf(Z, U, V, W_aug, nuep);
            W_hist{l_j + 1} = W_aug;
        end

        V(j, (L + 1):(L + n_new_lim)) = 0;
        trunc = exp(trunc - max(trunc));
        trunc = trunc/sum(trunc);
        p = rand;
        t = 0;
        for l_j=0:n_new_lim
            t = t+trunc(l_j+1);
            if p < t
                new_dishes = l_j;
                break;
            end;
        end;
        V(j, (L + 1):(L + new_dishes)) = 1; % exchangeable
        V = V(:, 1:(L + new_dishes)); % shink the zero columns
        W = W_hist{new_dishes + 1};
        L = L + new_dishes;
    end
end






















































