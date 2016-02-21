function output = mucondpdf(mur, murp, murm, a, b, c)
% Conditional probability Prob(mur | murp, murm) in Eq. (2.20)
% Referred to Eq. (19) in Xuan et. al [2015]
% Note: 1) mu's are TWO-DIM when paired and ONE-DIM when not!
%       2) mu's dimensions have to be identical
%       2) only need bi-variate pdf

% Arguments:
% mur -- mu_r
% murp -- mu_{r+1}
% murm -- mu_{r-1}
% a, b, c -- parameters of bi-variate beta distribution

% Todo:
%   - verify the consistency of -1
%   - verify murp and murm not zeros at same time
    
    % check if mu contains -1 that indicates non-pairing
    % if yes, the dimension should be all 1's for mur, murp, murm
    if ~any(mur == -1)
        if sum(murm) == 2 % starting mu (both 1)
            difmup = murp ./ mur;
            output = bibetapdf(difmup(1), difmup(2), a, b, c) / prod(mur);
        elseif sum(murp) == 0 % trailing my (both 0)
            difmum = mur ./ murm;
            output = bibetapdf(difmum(1), difmum(2), a, b, c) / prod(murm);
        else
            % both muru and murv exist and bibeta(a, b, c) is used
            difmum = mur ./ murm;
            difmup = murp ./ mur;
            output = bibetapdf(difmum(1), difmum(2), a, b, c) * ...
                bibetapdf(difmup(1), difmup(2), a, b, c) / ...
                prod(murm) / prod(mur);
        end
    elseif find(mur == -1) == 1
        % muru is null and beta(b, c) is used
        difmum = mur(2)/ murm(2);
        difmup = murp(2)/ mur(2);
        output = betapdf(difmum, b, c) * betapdf(difmup, b, c) / ... 
            murm(2) / mur(2);
    elseif find(mur == -1) == 2
        % murv is null and beta(a, c) is used
        difmum = mur(1)/ murm(1);
        difmup = murp(1)/ mur(1);
        output = betapdf(difmum, a, c) * betapdf(difmup, a, c) / ... 
            murm(1) / mur(1);
    end
end