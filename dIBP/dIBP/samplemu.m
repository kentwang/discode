function [mu_u, mu_v] = samplemu(U, V, K_plus, L_plus, mu_u, mu_v, a, b)
% Sample feature generating probabilities using MH algorithm

% INPUT
%   - U, V, mu_u, mu_v: current real feature matrices and probabilities
%   - a, b: controling parameters of dIBP

    I = size(U, 1);
    J = size(V, 1);

    % expan mu_u and mu_v when K_plus and L_plus are increased
    % or shink mu_u and mu_v when K_plus and L_plus are decreased

    K_u = length(mu_u);
    L_v = length(mu_v);

    if K_u > K_plus
        mu_u = mu_u(1:K_plus);
    elseif K_u < K_plus
        for k = (K_u+1):K_plus
            mu_u = [mu_u; mu_u(k-1)/3];
        end
    end

    if L_v > L_plus
        mu_v = mu_v(1:L_plus);
    elseif L_v < L_plus
        for l = (L_v+1):L_plus
            mu_v = [mu_v; mu_v(l-1)/3];
        end
    end


    M_low = min(K_plus, L_plus);
    M_high = max(K_plus, L_plus);

    % M-H sampling for paired mu. % M_low >= 2
    for r = 1:M_low
        mu_r = [mu_u(r), mu_v(r)];
        if r == 1
            mu_r_plus = [mu_u(r+1), mu_v(r+1)];

            mu_u_r_prop = truncbetarnd(a / K_plus, 1, mu_r_plus(1), 1);
            mu_v_r_prop = truncbetarnd(b / L_plus, 1,mu_r_plus(2), 1);
            mu_r_prop = [mu_u_r_prop, mu_v_r_prop];

            logl_curr = binologlik(U(:, r), mu_r(1)) + binologlik(V(:, r), mu_r(2)) + mu1_mu2logpdf(mu_r, mu_r_plus, a, b) - mu1_mu2proplogpdf(mu_r, mu_r_plus, a/K_plus, b/L_plus);
            logl_prop = binologlik(U(:, r), mu_r_prop(1)) + binologlik(V(:, r), mu_r_prop(2)) + mu1_mu2logpdf(mu_r_prop, mu_r_plus, a, b) - mu1_mu2proplogpdf(mu_r_prop, mu_r_plus, a/K_plus, b/L_plus);
        elseif r == M_low
            mu_r_minus = [mu_u(r-1), mu_v(r-1)];

            mu_u_r_prop = truncbetarnd(a / K_plus, 1, 0, mu_r_minus(1));
            mu_v_r_prop = truncbetarnd(b / L_plus, 1, 0, mu_r_minus(2));
            mu_r_prop = [mu_u_r_prop, mu_v_r_prop];

            logl_curr = binologlik(U(:, r), mu_r(1)) + binologlik(V(:, r), mu_r(2)) + mu2_mu1logpdf(mu_r_minus, mu_r, a, b) - mu2_mu1proplogpdf(mu_r_minus, mu_r, a/K_plus, b/L_plus);
            logl_prop = binologlik(U(:, r), mu_r_prop(1)) + binologlik(V(:, r), mu_r_prop(2)) + mu2_mu1logpdf(mu_r_minus, mu_r_prop, a, b) - mu2_mu1proplogpdf(mu_r_minus, mu_r_prop, a/K_plus, b/L_plus);
        else
            mu_r_plus = [mu_u(r+1), mu_v(r+1)];
            mu_r_minus = [mu_u(r-1), mu_v(r-1)];

            mu_u_r_prop = truncbetarnd(a / K_plus, 1, mu_r_plus(1), mu_r_minus(1));
            mu_v_r_prop = truncbetarnd(b / L_plus, 1, mu_r_plus(2), mu_r_minus(2));
            mu_r_prop = [mu_u_r_prop, mu_v_r_prop];

            logl_curr = binologlik(U(:, r), mu_r(1)) + binologlik(V(:, r), mu_r(2)) + mu_fullcondlogpdf(mu_r_minus, mu_r, mu_r_plus, a, b) - truncbetalogpdf(mu_r(1), a/K_plus, 1, mu_r_plus(1), mu_r_minus(1)) - truncbetalogpdf(mu_r(2), b/L_plus, 1, mu_r_plus(2), mu_r_minus(2));
            logl_prop = binologlik(U(:, r), mu_r_prop(1)) + binologlik(V(:, r), mu_r_prop(2)) + mu_fullcondlogpdf(mu_r_minus, mu_r, mu_r_plus, a, b) - truncbetalogpdf(mu_r_prop(1), a/K_plus, 1, mu_r_plus(1), mu_r_minus(1)) - truncbetalogpdf(mu_r_prop(2), b/L_plus, 1, mu_r_plus(2), mu_r_minus(2));
        end

        % M-H updating
        mu_acc = exp(min(0, logl_prop - logl_curr));
        % disp(['accept rate: ', num2str(mu_acc)]);

        if rand < mu_acc
            % disp('updating...')
            mu_u(r) = mu_r_prop(1);
            mu_v(r) = mu_r_prop(2);
        end
    end

    % M-H sampling for possible unpaired mu using univariate
    if K_plus < M_high % update mu_v
        for r = (M_low+1):M_high
            if r == M_high
                mu_v_r_minus = mu_v(r-1);
                mu_v_r_prop = truncbetarnd(b/L_plus, 1, 0, mu_v_r_minus);

                logl_curr = binologlik(V(:, r), mu_v(r)) + mu2_mu1logpdf_univ(mu_v_r_minus, mu_v(r), b) - mu2_mu1proplogpdf_univ(mu_v_r_minus, mu_v(r), b/L_plus);

                logl_prop = binologlik(V(:, r), mu_v_r_prop) + mu2_mu1logpdf_univ(mu_v_r_minus, mu_v_r_prop, b) - mu2_mu1proplogpdf_univ(mu_v_r_minus, mu_v_r_prop, b/L_plus);
            else
                mu_v_r_minus = mu_v(r-1);
                mu_v_r_plus = mu_v(r+1);
                mu_v_r_prop = truncbetarnd(b/L_plus, 1, mu_v_r_plus, mu_v_r_minus);

                logl_curr = binologlik(V(:, r), mu_v(r)) + mu_fullcondlogpdf_univ(mu_v_r_minus, mu_v(r), mu_v_r_plus, b) - truncbetalogpdf(mu_v(r), b/L_plus, 1, mu_v_r_plus, mu_v_r_minus);

                logl_prop = binologlik(V(:, r), mu_v_r_prop) + mu_fullcondlogpdf_univ(mu_v_r_minus, mu_v_r_prop, mu_v_r_plus, b) - truncbetalogpdf(mu_v_r_prop, b/L_plus, 1, mu_v_r_plus, mu_v_r_minus);
            end

            % update unmatched mu_v
            mu_v_acc = exp(min(0, logl_prop - logl_curr));

            if rand < mu_v_acc
                mu_v(r) = mu_v_r_prop;
            end   
        end
    elseif L_plus < M_high % update mu_u. THIS HAS NOT BEEN TESTES YET.
        for r = (M_low+1):M_high
            if r == M_high
                mu_u_r_minus = mu_u(r-1);
                mu_u_r_prop = truncbetarnd(a/K_plus, 1, 0, mu_u_r_minus);
                % disp(mu_u_r_prop);

                logl_curr = binologlik(U(:, r), mu_u(r)) + mu2_mu1logpdf_univ(mu_u_r_minus, mu_u(r), a) - mu2_mu1proplogpdf_univ(mu_u_r_minus, mu_u(r), a/K_plus);
                % disp(logl_curr);

                logl_prop = binologlik(U(:, r), mu_u_r_prop) + mu2_mu1logpdf_univ(mu_u_r_minus, mu_u_r_prop, a) - mu2_mu1proplogpdf_univ(mu_u_r_minus, mu_u_r_prop, a/K_plus);
                % disp(logl_prop - logl_curr);

                % mu_u_acc = exp(min(0, logl_prop - logl_curr));
                % disp(['accept rate: ', num2str(mu_u_acc)]);

            else
                mu_u_r_minus = mu_u(r-1);
                mu_u_r_plus = mu_u(r+1);
                mu_u_r_prop = truncbetarnd(a/K_plus, 1, mu_u_r_plus, mu_u_r_minus);

                logl_curr = binologlik(U(:, r), mu_u(r)) + mu_fullcondlogpdf_univ(mu_u_r_minus, mu_u(r), mu_u_r_plus, b) - truncbetalogpdf(mu_u(r), a/K_plus, 1, mu_u_r_plus, mu_u_r_minus);

                logl_prop = binologlik(U(:, r), mu_u_r_prop) + mu_fullcondlogpdf_univ(mu_u_r_minus, mu_u_r_prop, mu_u_r_plus, a) - truncbetalogpdf(mu_u_r_prop, a/K_plus, 1, mu_u_r_plus, mu_u_r_minus);

                % mu_u_acc = exp(min(0, logl_prop - logl_curr));
                % disp(['accept rate: ', num2str(mu_u_acc)]);
            end

            % update unmatched mu_u
            mu_u_acc = exp(min(0, logl_prop - logl_curr));
            % disp(['accept rate: ', num2str(mu_u_acc)]);

            if rand < mu_u_acc
                mu_u(r) = mu_u_r_prop;
                % disp('updating...')
            end 
        end
    else
        continue;
    end
end