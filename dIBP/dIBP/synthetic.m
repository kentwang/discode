function [U, V, Z, X] = synthetic(W, I, J, K, L, nuep)
% Simulate ground truth data
% Output:
%   - U, V,
%   - Z, X
    
    U = featrnd(I, K, 0.7);
    V = featrnd(J, L, 0.8);
    
    Eps = normrnd(nuep, 1, I, J);
    Z = U * W * V' + Eps;
    X = (Z > 0);
end