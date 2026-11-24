function [P, dP] = legendre_basis(k, s)
%LEGENDRE_BASIS  Legendre basis of degree up to k on [0, 1].
%
%   [P, dP] = legendre_basis(k, s) evaluates the k+1 functions and their
%   derivatives at the column vector s.
    arguments
        k double
        s (:, 1) double
    end
    d1 = k + 1;
    [P, dP] = legendre_eval(1:d1, s);
end