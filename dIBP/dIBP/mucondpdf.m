function output = mucondpdf(mur, murp, murm, a, b, c)
% Conditional probability Prob(mur | murp, murm) in Eq. (2.20)
% Referred to Eq. (19) in Xuan et. al [2015]
% Note: 1) all mu's are TWO-DIM!
%       2) only need bi-variate pdf

% Arguments:
% mur -- mu_r
% murp -- mu_{r+1}
% murm -- mu_{r-1}
% a, b, c -- parameters of bi-variate beta distribution
    
    difmum = mur ./ murm;
    difmup = murp ./ mur;
    output = bibetapdf(difmum(1), difmum(2), a, b, c) * ...
        bibetapdf(difmup(1), difmup(2), a, b, c) / ...
        prod(murm) / prod(mur);
end