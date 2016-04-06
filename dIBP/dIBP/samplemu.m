function [mu_u, mu_v] = samplemu(U, V, mu_u, mu_v, a, b)
% Sample feature generating probabilities using MH algorithm

% INPUT
%   - U, V, mu_u, mu_v: current real feature matrices and probabilities
%   - a, b: controling parameters of dIBP
    
% TODO
%   - mu_u is matched with U

    [I, K] = size(U);
    [J, L] = size(V);
    M_low = min(K, L)
    M_high = max(K, L);

    % M-H sampling for paired mu. % M_low >= 2
    for r = 1:M_low
        mu_r = [mu_u(r), mu_v(r)];
        if r == 1
            mu_r_plus = [mu_u(r+1), mu_v(r+1)];

            mu_u_r_prop = truncbetarnd(a / K, 1, mu_u(r+1), 1);
            mu_v_r_prop = truncbetarnd(b / L, 1, mu_v(r+1), 1);
            mu_r_prop = [mu_u_r_prop, mu_v_r_prop];

            logl_curr = binologlik(U(:, r), mu_r(1)) + binologlik(V(:, r), mu_r(2)) + mu1_mu2logpdf(mu_r, mu_r_plus, a, b) - mu1_mu2proplogpdf(mu_r, mu_r_plus, a/K, b/L);
            logl_prop = binologlik(U(:, r), mu_r_prop(1)) + binologlik(V(:, r), mu_r_prop(2)) + mu1_mu2logpdf(mu_r_prop, mu_r_plus, a, b) - mu1_mu2proplogpdf(mu_r_prop, mu_r_plus, a/K, b/L);
        elseif r == M_low
            mu_r_minus = [mu_u(r-1), mu_v(r-1)];

            mu_u_r_prop = truncbetarnd(a / K, 1, 0, mu_u(r-1));
            mu_v_r_prop = truncbetarnd(b / L, 1, 0, mu_v(r-1));
            mu_r_prop = [mu_u_r_prop, mu_v_r_prop];

            logl_curr = binologlik(U(:, r), mu_r(1)) + binologlik(V(:, r), mu_r(2)) + mu2_mu1logpdf(mu_r_minus, mu_r, a, b) - mu2_mu1proplogpdf(mu_r_minus, mu_r, a/K, b/L);
            logl_prop = binologlik(U(:, r), mu_r_prop(1)) + binologlik(V(:, r), mu_r_prop(2)) + mu2_mu1logpdf(mu_r_minus, mu_r_prop, a, b) - mu2_mu1proplogpdf(mu_r_minus, mu_r_prop, a/K, b/L);

        else
        end

        % M-H updating
        mu_acc = exp(min(0, logl_prop - logl_curr));
        if rand < mu_acc
            mu(r, :) = mu_r_prop;
        end        
    end
    
    % append zerors for U or V
    % if K < M
    %     U = [U, zeros(I, M - K)];
    % else
    %     V = [V, zeros(J, M - J)];
    % end

    % % append zeros for mu_u and mu_v based on size of U/V
    % if length(mu_u) < M
    %     for i = (length(mu_u) + 1):M
    %         mu_u = [mu_u; mu_u(i - 1) / 3];
    %     end
    % end

    % if length(mu_v) < M
    %     for i = (length(mu_v) + 1):M
    %         mu_v = [mu_v; mu_v(i - 1) / 3];
    %     end
    % end

    % mu = [mu_u, mu_v];
    
    % MH sampling
    for r = 1:M
        % input of each sampling        
        if r == 1 % propose for the first mu
            lower_u = mu_u(r + 1);
            lower_v = mu_v(r + 1);
            upper_u = min(mu_u(r) + 0.1, 1);
            upper_v = min(mu_v(r) + 0.1, 1);
        elseif r == M % propose for the last mu
            lower_u = max(mu_u(r) - 0.1, 0);
            lower_v = max(mu_v(r) - 0.1, 0);
            upper_u = mu_u(r - 1);
            upper_v = mu_v(r - 1);
        else % update others
            lower_u = mu_u(r + 1);
            lower_v = mu_v(r + 1);
            upper_u = mu_u(r - 1);
            upper_v = mu_v(r - 1);
        end
        
        mu_u_r_prop = truncbetarnd(a / K, 1, lower_u, upper_u);
        mu_v_r_prop = truncbetarnd(b / L, 1, lower_v, upper_v);
        mu_r_prop = [mu_u_r_prop, mu_v_r_prop];
        
        if r == 1 % starting mu accept ratio
            accept_r = muacceptratio(U(:, r), V(:, r), mu(r, :), ...
                mu_r_prop, mu(r + 1, :), [1, 1], a, b, 1, K, L);
        elseif r == M % trailing mu accept ratio
            accept_r = muacceptratio(U(:, r), V(:, r), mu(r, :), ...
                mu_r_prop, [0, 0], mu(r - 1, :), a, b, 1, K, L);
        else
            accept_r = muacceptratio(U(:, r), V(:, r), mu(r, :), ...
                mu_r_prop, mu(r + 1, :), mu(r - 1, :), a, b, 1, K, L);
        end
        
        if rand < accept_r
            mu(r, :) = mu_r_prop;
        end
    end
    
    mu_u = mu(1:K, 1); % truncated to real probs
    mu_v = mu(1:L, 2);
end



