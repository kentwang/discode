function [U, V, K_plus, L_plus] = sampleUV(Z, U, V, I, J, K_plus, L_plus, nuep, a, b, sigma_w)
% Update U and V (simultaneously K and L)
% Augment W for new features of U and V

% Todo
%   - Add K_plus and L_plus as input
%   - Fix the updating of K_plus and L_plus
%   - Integrate out W in the very beginning

% Test
% U = U_true; V = V_true; W = W_true; Z = Z_true;
    
    % update U (and W if any new row features from U)
    % [I, K] = size(U);
    for i = 1:I
        % printf('U iter %d/%d\n', i, I);
        for k = 1:K_plus
            if k > K_plus
                break;
            end
            
            % shrink singular column
            % Note: if all zeros except U(i, k), then p(u_{ik}=1|) = 0!
            % Todo: - what if k = K  and singular?
            if U(i, k) > 0
                if sum(U(:, k)) - U(i, k) == 0
                    U(i, k) = 0;
                    U(:, k:(K_plus-1)) = U(:, (k + 1):K_plus);
                    K_plus = K_plus - 1;
                end
            end
                    
            % sample U
            U(i, k) = 1;
            p(1) = Z_UVlogpdf(Z, U(:,1:K_plus), V(:,1:L_plus), sigma_w, nuep, I, J, K_plus, L_plus) + log(sum(U(:, k) - 1)) - log(I);
            U(i, k) = 0;
            p(2) = Z_UVlogpdf(Z, U(:,1:K_plus), V(:,1:L_plus), sigma_w, nuep, I, J, K_plus, L_plus) + log(I - sum(U(:, k))) - log(I);
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
        % W_hist = {};
        alpha_I = a / I; % alpha/I
        
        for k_i = 0:n_new_lim
            U(i, (K_plus+1):(K_plus+k_i)) = 1; % test here U with new feature
            % W_aug = W_expand(k_i, 2, Z, U, V, W, nuep, sigma_w); % append rows for new features U. TEST GOOD!
            trunc(k_i+1) = k_i*log(alpha_I) - alpha_I - log(factorial(k_i)) + Z_UVlogpdf(Z, U(:,1:K_plus+k_i), V, sigma_w, nuep, I, J, K_plus+k_i, L_plus);
            % W_hist{k_i + 1} = W_aug;
        end
        
        U(i, (K_plus+1):(K_plus+n_new_lim)) = 0;
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
        U(i, (K_plus+1):(K_plus+new_dishes)) = 1; % exchangeable
        % U = U(:, 1:(K + new_dishes)); % shink the zero columns. NO NEED!
        % W = W_hist{new_dishes + 1};
        K_plus = K_plus + new_dishes;
    end

    % plot the updated U matrix
    % clf; imagesc(U);

    % update V (and W if any new column features from V)
    % [J, L] = size(V);
    for j = 1:J
        % printf('V iter %d/%d\n', j, J);
        for l = 1:L_plus
            if l > L_plus
                break;
            end

            if V(j, l) > 0
                if sum(V(:, l)) - V(j, l) == 0
                    V(j, l) = 0;
                    V(:, l:(L_plus-1)) = V(:, (l + 1):L_plus);
                    L_plus = L_plus - 1;
                end
            end

            % sample V
            V(j, l) = 1;
            p(1) = Z_UVlogpdf(Z, U(:,1:K_plus), V(:,1:L_plus), sigma_w, nuep, I, J, K_plus, L_plus) + log(sum(V(:, l) - 1)) - log(J);
            V(j, l) = 0;
            p(2) = Z_UVlogpdf(Z, U(:,1:K_plus), V(:,1:L_plus), sigma_w, nuep, I, J, K_plus, L_plus) + log(J - sum(V(:, l))) - log(J);
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
            V(j, (L_plus+1):(L_plus+l_j)) = 1; % test here V with new feature
            % W_aug = W_expand(l_j, 1, Z, U, V, W, nuep, sigma_w);
            trunc(l_j+1) = l_j * log(beta_J) - beta_J - log(factorial(l_j)) + Z_UVlogpdf(Z, U(:,1:K_plus), V(:,1:L_plus+l_j), sigma_w, nuep, I, J, K_plus, L_plus+l_j);
            % W_hist{l_j + 1} = W_aug;
        end

        V(j, (L_plus+1):(L_plus+n_new_lim)) = 0;
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
        V(j, (L_plus+1):(L_plus+new_dishes)) = 1; % exchangeable
        % V = V(:, 1:(L + new_dishes)); % shink the zero columns
        % W = W_hist{new_dishes + 1};
        L_plus = L_plus + new_dishes;
    end
end