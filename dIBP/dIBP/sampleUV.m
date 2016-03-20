function U = sampleUV(U, V, W, nuep, a)
% Update U and V (simultaneously K and L)

Test
U = U_true; V = V_true;
W = W_true; Z = Z_true;

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
            U_1 = U; U(i, k) = 1;
            p(1) = Z_UVWlogpdf(U_1, V, W, nuep) + log(sum(U(:, k) - 1)) - ...
                log(I);
            U_0 = U; U(i, k) = 0;
            p(2) = Z_UVWlogpdf(U_0, V, W, nuep) + log(I - sum(U(:, k))) - ...
                log(I);
            p = exp(p - max(p));
            
            if rand < p(1)/(p(1)+p(2))
                U(i,k) = 1;
            else
                U(i,k) = 0;
            end
        end

        % sample new features
        n_new_lim = 4;
        trunc = zeros(1, n_new_lim + 1);
        W_hist = {};
        alpha_I = a / I; % alpha/I
        
        for k_i = 0:n_new_lim
            U(i, (K + 1):(K + k_i)) = 1; % test here U with new feature
            W_aug = W_expand(k_i, 2, U, V, W, nuep, sigma_w); % append rows for new features U. TEST GOOD!
            trunc(k_i + 1) = k_i * log(alpha_I) - alpha_I - ...
                log(factorial(k_i)) + Z_UVWlogpdf(U, V, W_aug, nuep);
            W_hist{k_i + 1} = W_aug;
        end
        
        U(i, (K + 1):(K + n_new_lim)) = 0;
        trunc = exp(trunc - max(trunc));
        trunc = trunc/sum(trunc);
        p = rand;
        t = 0;
        for k_i=0:4
            t = t+trunc(k_i+1);
            if p < t
                new_dishes = k_i;
                break;
            end;
        end;
        U(i, (K + 1):(K + new_dishes)) = 1; % exchangeable
        W = W_hist{new_dishes + 1}
        K = K + new_dishes;
    end

end