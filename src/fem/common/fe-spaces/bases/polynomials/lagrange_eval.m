function [P, dPx, dPy] = lagrange_eval(j, r, s)
%-------------------------------------------------------------------------
%Description: Evaluates the Lagrange basis functions and their first
%           derivatives at the natural coordinates (r, s). Only j = 1 is
%           implemented: the three barycentric functions r, s, 1 - r - s.
%Input:     - j: scalar, polynomial degree
%           - r, s: column vectors of natural coordinates, same length
%Output:    - P: basis function values, one column per function
%           - dPx, dPy: first partial derivatives w.r.t. r and s
%Notes:     j other than 1, and r and s of different length, raise an
%           error.
%-------------------------------------------------------------------------
    arguments
        j (1, 1) double
        r (:, 1) double
        s (:, 1) double
    end

    if j ~= 1
        error("lagrange_eval:degree", "Only degree 1 is implemented.")
    end

    if length(r) ~= length(s)
        error("Natural coordinates do not match.")
    end

    % barycentric coordinates
    L = [r, s, 1 - r - s];
    
    [P, dPx, dPy] = deal(zeros(length(r), length(j)));

    if j == 1
        P = L;
        dPx = L.^0 .* [1, 0, -1];
        dPy = L.^0 .* [0, 1, -1];
    end
    
end