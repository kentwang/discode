function output = Z_UVWlogpdf(Z, U, V, W, nuep)
% log conditional pdf of Z on U, V, W, nuep

% Test
% U = U_true; V = V_true;
% W = W_true; Z = Z_true;

    output = sum(sum(log(normpdf(Z, U * W * V' + nuep, 1))));
end