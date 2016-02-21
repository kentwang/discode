function output = samplemu(U, V, mu_u, mu_v, a, b)
% Sample feature generating probabilities using MH algorithm

% INPUT
%   - U, V, mu_u, mu_v: current real feature matrices and probabilityes
%   - a, b: controling parameters of dIBP
    
% TODO
%   - Warning unmatched dimensions
%   - recheck 0 and 1 as the bounds. Maybe a noise +/- is OK

    K = length(mu_u);
    L = length(mu_v);
    M = max(K, L);
    
    % append zerors
    if(K < M) 
        U = [U, zeros(I, M - K)];
        mu_u = [mu_u; zeros(M - K, 1)];
    else
        V = [V, zeros(J, M - J)];
        mu_v = [mu_v; zeros(M - L, 1)];
    end
    
    % MH sampling
    for r = 1:M
        % input of each sampling        
        if r == 1 % propose for the first mu
            lower_u = mu_u(r + 1);
            lower_v = mu_v(r + 1);
            upper_u = 1;
            upper_v = 1;
        elseif r == M % propose for the last mu
            lower_u = 0;
            lower_v = 0;
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
        
        accept_r = min(1, UV_mupdf(U(:, r), V(:, r), mu_u_r_prop, ...
            mu_v_r_prop) / UV_mupdf(U(:, r), V(:, r), mu_u(r), ...
            mu_v(r)) * mucondpdf(mur_test, murp_test, murm_test, ...
            a, b, c) / mucondpdf(mur, murp, murm, a, b, c) * ...
            mucondproppdf(mur, murp, murm, a, b, c, K, L) / ...
            mucondproppdf(mur_test, murp_test, murm_test, a, b, c, K, L));
    end
end



