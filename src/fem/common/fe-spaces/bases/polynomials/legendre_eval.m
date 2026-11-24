function [P, dP] = legendre_eval(j_list, s)
%-------------------------------------------------------------------------
%Description: Evaluates the Legendre basis function of index j_list,
%           and its first derivative, at the natural coordinate s in [0,1]
%Input:     - j_list: basis function index (>= 1)
%           - s: natural coordinate
%Output:    - P: basis function values
%           - dP: first derivative w.r.t. s
%-------------------------------------------------------------------------
arguments
    j_list (1, :) double
    s (:, :, :) double
end
s = reshape(s(:), [], size(s, 3));

t = 2*s - 1;
dtds = 2;

[P, dP] = deal(zeros(size(s, 1), length(j_list), size(s, 2)));
for i = 1:length(j_list); j = j_list(i);
    n = j - 1;
    [p1, dpdt] = jacobi_eval(n,0,0,t);
    p = p1;
    dp = dpdt * dtds;

    np = sqrt(1/(2*n+1));
    P(:, i, :) = p/np; 
    dP(:, i, :) = dp/np;
end
end