function output = samplemu(U, V, mu_u, mu_v)
% Sample feature generating probabilities using MH algorithm

% INPUT
%   - U, V, mu_u, mu_v: current feature matrices and probabilityes

% NOTE
%   - Dimensions of U (V) and mu_u (mu_v) have to be matched
    
% TODO
%   - Warning unmatched dimensions

    K = length(mu_u);
    L = length(mu_v);
    Mm = min(K, L);
    Mp = max(K, L);
    
    
end